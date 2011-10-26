/*
 *  Copyright (C) 2010 by Tim Massingham
 *  tim.massingham@ebi.ac.uk
 *
 *  This file is part of the AYB base-calling software.
 *
 *  AYB is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AYB is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AYB.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <omp.h>
#include <assert.h>
#include <math.h>
#include <strings.h>
#include <stdlib.h>
#include "matrix.h"
#include "nuc.h"
#include "process_intensities.h"
#include "ayb.h"
#include "call_bases.h"
#include "options.h"
#include "lapack.h"

#include CALIBRATION_FILE

const static real_t polyQual = 38.0;

real_t logadd(const real_t a, const real_t b){
	return (a>b)?
	       a + log1p(exp(b-a)) :
	       b + log1p(exp(a-b));
}

/**  Subcalculation for call_bases
 *  om is pointer to first entry of i^th block on the diagonal of Omega. Lda of Omega is 4*ncycle
 *  p is intensities of i^th cycle
 *  b is putative base
 *  Returns -2*I_b^t om p + lambda * I_b^t om I^b, where I_b is indicator vector of base b
 */
real_t baselike (const real_t * restrict p, const real_t lambda, const NUC b, const real_t * om, const int ncycle){
       if(NULL==p || NULL==om || !isfinite(lambda) || NUC_AMBIG==b){ return NAN;}
       const int lda = 4*ncycle;
       real_t like = lambda*om[b*lda+b];
       for ( int i=0 ; i<4 ; i++){
               like -= 2.0* p[i]*om[b*lda+i];
       }
       return like;
}

/**  Subcalculation for call_bases
 *  om is pointer to first entry of i^th block on the lower diagonal of Omega. Lda of Omega is 4*ncycle
 *  p is intensities of i^th cycle
 *  a is putative base for cycle i
 *  b is putative base for cycle i+1
 *  Returns lambda I_b^t om I_a - p_{i+1}^t om p_{i} - p_{i}^t om p_{i+1}^t
 */
real_t crosslike(const real_t * restrict p, const real_t lambda, const NUC a, const NUC b, const real_t * om, const int ncycle){
	if(NULL==p || NULL==om || !isfinite(lambda) || NUC_AMBIG==b || NUC_AMBIG==b){ return NAN;}
        const int lda = 4*ncycle;
        real_t like = lambda*om[a*lda+b];
        for ( int i=0 ; i<4 ; i++){
                like -= p[i]*om[i*lda+b];
                like -= p[4+i]*om[a*lda+i];
        }
        return like;
}

/**  Call bases from processed intensities
 *  Uses a dynamic programming algorithm (Viterbi) to find the bases that minimises the 
 *  squared error: (p-lambda*s)^t Omega (p-lambda*s).
 *               = p^t Om p + lambda ( -2 * s^t Omega p + lambda s^t Omega s )
 *                                   |------------------ A -------------------|
 *  Sufficient to minimise A
 *  Omega is block-tridiagonal, that is composed from 4x4 blocks, e.g. Om_{11} Om_{21}^t 0        
 *                                                                     Om_{21} Om_{22}   Om_{32}^t
 *                                                                     0       Om_{32}   Om_{33}  
 */  
real_t call_base( const MAT p, const real_t lambda, const MAT omega, NUC * base){
	if(NULL==base || NULL==p || NULL==omega){ return NAN; }
        const int ncycle = p->ncol;
        const int lda = omega->ncol;
       
        // array contains accumulation information and trace-back directions
        real_t array[4*ncycle];
        // Initialise. First elements of array contain contribution from first cycle.       
        for ( int b=0 ; b<NBASE ; b++){
                array[b] = baselike(p->x,lambda,b,omega->x,ncycle);
        }

        // Forwards piece of algorithm. For each cycle, for each base, find the base in the
        // previous cycle that minimises 
        for ( int cy=1; cy<ncycle ; cy++){
                NUC precall[NBASE];
                for ( int b=0 ; b<NBASE ; b++){
                	// Contribution from calling b at cycle cy, independent from other cycles
                        real_t diag = baselike(p->x+cy*NBASE,lambda,b,omega->x+cy*NBASE*lda+cy*NBASE,ncycle);
                        // Find call prev at previous cycle that minimises cost of calling b and prev
                        real_t minstat = HUGE_VAL;
                        int minidx = 0;
                        real_t stat[NBASE];
                        for ( int prev=0 ; prev<NBASE ; prev++){
                                stat[prev] = diag +                       // Cost of calling b at cycle cy
                                	     array[(cy-1)*NBASE+prev] +   // Cost of calling prev at cycle (cy-1)
                                	     // Adjustment for calling both.
                                	     2.0*crosslike(p->x+(cy-1)*NBASE,lambda,prev,b,omega->x+(cy-1)*NBASE*lda+cy*NBASE,ncycle);
                                // Keep track of best previous call
                                if(stat[prev]<minstat){
                                        minidx = prev;
                                        minstat = stat[prev];
                                }
                        }
			// Save best call for previous cycle given call b at this cycle.
                        array[cy*NBASE+b] = minstat;
                        precall[b] = minidx;
                }
                // No longer need previous entries of array. Use memory to save trace-back 
                // information. 
                for ( int b=0 ; b<NBASE ; b++){
                        array[(cy-1)*NBASE+b] = precall[b];
                }
        }
        
        // Backwards piece of algorithm.
        // Find best call on last cycle
        real_t minstat = HUGE_VAL;
        int minidx = 0;
        for ( int b=0 ; b<NBASE ; b++){
                if(array[(ncycle-1)*NBASE+b]<minstat){ minidx = b; minstat = array[(ncycle-1)*NBASE+b];}
        }
	// Trace calls backward using previously stored information
        base[ncycle-1] = minidx;
        for ( int cy=(ncycle-2) ; cy>=0 ; cy--){
                base[cy] = (NUC)array[cy*NBASE+base[cy+1]];
        }
	return xMy(p->x,omega,p->x) + lambda * minstat;
}

/** Posterior probabilities via fwds/bwds
 *  Repeats much of call_bases
 */
void call_qualities_post(const MAT p, const real_t lambda, const MAT omega, const real_t effDF, NUC * base, real_t * qual){
	if(NULL==base || NULL==p || NULL==omega){ return; }
        const int ncycle = p->ncol;
        const int lda = omega->ncol;
	const real_t polyErr = -expm1(-polyQual/10.0 * log(10.0));

        // arrays contain accumulation information
        real_t farray[4*ncycle], barray[4*ncycle];
        // Calculate costs
	real_t basecost[4*ncycle], crosscost[16*ncycle];
	for ( int b=0 ; b<NBASE ; b++){
		basecost[b] = baselike(p->x,lambda,b,omega->x,ncycle);
	}
	for ( int cy=1 ; cy<ncycle ; cy++){
        	for ( int b=0 ; b<NBASE ; b++){
                	basecost[4*cy+b] = baselike(p->x+cy*NBASE,lambda,b,omega->x+cy*NBASE*lda+cy*NBASE,ncycle);
			for ( int b2=0 ; b2<NBASE ; b2++){
				crosscost[16*cy+4*b+b2] = crosslike(p->x+(cy-1)*NBASE,lambda,b,b2,omega->x+(cy-1)*NBASE*lda+cy*NBASE,ncycle);
			}
		}
        }

        // Forwards piece of algorithm.
	// Initalise
	for ( int b=0 ; b<4 ; b++){ farray[b] = 0.0; }
	// Iteration
        for ( int cy=1; cy<ncycle ; cy++){
                for ( int b=0 ; b<NBASE ; b++){
			farray[4*cy+b] = HUGE_VAL;
			// Find minimum for past calls give current call
			for ( int prev=0 ; prev<4 ; prev++){
				real_t callLS = basecost[4*(cy-1)+prev] + farray[4*(cy-1)+prev] + 2.0*crosscost[16*cy+4*prev+b];
				if(callLS<farray[4*cy+b]){ farray[4*cy+b] = callLS; }
			}
                }
        }
	real_t fwdslike = HUGE_VAL;
	for ( int b=0 ; b<4 ; b++){
		if(basecost[(ncycle-1)*4+b]+farray[(ncycle-1)*4+b]<fwdslike){
			fwdslike = basecost[(ncycle-1)*4+b]+farray[(ncycle-1)*4+b];
		}
	}


        // Backwards piece of algorithm.
        // Initialise
	for ( int b=0 ; b<4 ; b++){ barray[(ncycle-1)*4+b] = 0.0; }
	// Iteration
	for ( int cy=(ncycle-2) ; cy>=0 ; cy--){
		for ( int b=0 ; b<NBASE ; b++){
			barray[4*cy+b] = HUGE_VAL;
			// Find minimum for future calls given current call
			for ( int nxt=0 ; nxt<4 ; nxt++){
				real_t callLS = basecost[4*(cy+1)+nxt] + barray[4*(cy+1)+nxt] + 2.0*crosscost[16*(cy+1)+4*b+nxt];
				if(callLS<barray[4*cy+b]){ barray[4*cy+b] = callLS; }
			}
		}
	}
	real_t bkdslike = HUGE_VAL;
	for ( int b=0 ; b<4 ; b++){
		if(basecost[b]+barray[b]<bkdslike){
			bkdslike = basecost[b]+barray[b];
		}
	}
	//fprintf(stdout,"%e %e\n",fwdslike,bkdslike);

	real_t pOp = xMy(p->x,omega,p->x);
	//fprintf(stdout,"%e\n",pOp+lambda*(farray[base[0]] + barray[base[0]] + basecost[base[0]]));
	for ( int cy=0 ; cy<ncycle ; cy++){
		int b = base[cy];
		real_t xmax = farray[4*cy+b] + barray[4*cy+b] + basecost[4*cy+b];
		xmax = pOp + lambda * xmax;
		real_t sum = 0.0;
		for ( int b=0 ; b<4 ; b++){
			real_t stat = farray[4*cy+b] + barray[4*cy+b] + basecost[4*cy+b];
			stat *= lambda;
			stat += pOp;
			sum += exp(-0.5*(1.0+effDF)*(log1p(stat)-log1p(xmax)));
		}
		real_t prob = 1.0 / sum;
		prob *= polyErr;
		qual[cy] = qual_from_prob(prob);
	}
}




real_t adjust_quality(const real_t qual, const NUC prior, const NUC base, const NUC next){

   real_t adj = calibration_intercept + calibration_scale * qual + calibration_baseprior_adj[prior*NBASE+base] + calibration_basenext_adj[next*NBASE+base] + calibration_priorbasenext_adj[(next*NBASE+prior)*NBASE+base];
   return adj;
}

real_t adjust_first_quality(const real_t qual, const NUC base, const NUC next){
   return calibration_intercept + calibration_scale * qual + calibration_basenext_adj[next*NBASE+base];
}

real_t adjust_last_quality(const real_t qual, const NUC prior, const NUC base){
   return calibration_intercept + calibration_scale * qual + calibration_baseprior_adj[prior*NBASE+base];
}


/* Calculates adjusted -loglikelihood from -loglikelihood
 * -log(mu+exp(-like))
 */
static inline real_t adjust_loglikelihood(const real_t negloglike, const real_t mu){
   const real_t prob = exp(-negloglike);
   return (prob>mu)?(negloglike-log1p(-mu/prob)):-(log(mu)-log1p(-prob/mu));
}

/*
 * Calculate likelihoods from processed intensities
 *  p       Processed intensities for given cycle
 *  lambda  Cluster brightness
 *  omega   Cycle specific inverse covariance matrix
 *  like    Vector of 4 reals for writing likelihoods to.
 *  Likelihoods are stored as -log(likelihood) == 0.5*least-squares statistic
 * REFACTOR: repeats much of call_base
 */
void call_likelihoods ( const real_t * restrict p, const real_t lambda, const MAT omega, real_t * like){
    assert(NULL!=p);
    assert(NULL!=omega);
    assert(NULL!=like);
    assert(NBASE==omega->nrow);
    assert(NBASE==omega->ncol);
    extern AYBOPT aybopt;

    if(0==lambda){
        for( int i=0 ; i<NBASE ; i++){
	   like[i] = (aybopt.mu)?-log(aybopt.mu):HUGE_VAL;
	}
	return;
    }

    const real_t K = xMy(p,omega,p);
    for ( int i=0 ; i<NBASE ; i++){
        like[i] = lambda * omega->x[i*NBASE+i];
        for ( int j=0 ; j<NBASE ; j++){
            like[i] -= 2.0 * p[j] * omega->x[i*NBASE+j];
        }
        like[i] = K + like[i]*lambda;
	like[i] /= 2.0;
    }

    if ( aybopt.mu ){
       for ( int i=0 ; i<NBASE ; i++){
	  like[i] = adjust_loglikelihood(like[i],aybopt.mu);
       }
    }
}




    
/* Maximum of n real_ts, should be moved into utility.h library */
static inline int max_real_t(const real_t * restrict p, const uint32_t n){
    validate(NULL!=p,-1);
    real_t m = p[0];
    int idx = 0;
    for ( uint32_t i=1 ; i<n ; i++){
        if(p[i]>m){ m=p[i]; idx=i;}
    }
    return idx;
}


/* 
 * Call base from processed intensities using maximum intensity
 */
NUC call_base_simple( const real_t * restrict p){
    return max_real_t(p,NBASE);
}


/*
 * Routines to calculate covariance of errors, using the "fast approach".
 * Working with processed intensities.
 */

/* Accumulate required variances (inner summation of variance calculation.
 * p        Matrix of processed intensities
 * lambda   Brightness of cluster
 * V        Array of covariance matrices into which acculation occurs.
 * Note: If V is NULL, the require memory is allocated
 */
MAT accumulate_covariance( const real_t we, const MAT p, const real_t lambda, const NUC * base, MAT V){
    validate(NULL!=p,NULL);
    validate(NBASE==p->nrow,NULL);
    validate(lambda>=0.,NULL);
    const uint32_t ncycle = p->ncol;
    const int lda = 4*ncycle;
    // Allocate memory for V if necessary
    if ( NULL==V ){
        V = new_MAT(lda,lda);
        if(NULL==V){ goto cleanup;}
    }
    // Perform accululation. V += R R^t
    // Note: R = P - \lambda I_b, where I_b is unit vector with b'th elt = 1
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        p->x[cy*NBASE+base[cy]] -= lambda;
    }
    const real_t alpha = 1.0;
    syr(LAPACK_LOWER,&lda,&alpha,p->x,LAPACK_UNIT,V->x,&lda);

    return V;
    
cleanup:
    free_MAT(V);
    return NULL;
}


/*
 *  Calculate covariance of (processed residuals) 
 */
MAT calculate_covariance( const AYB ayb){
    const uint32_t ncluster = ayb->ncluster;
    const uint32_t ncycle = ayb->ncycle;
    
    //MAT invAt = invert(ayb->At); //transpose(ayb->A);
    struct structLU AtLU = LUdecomposition(ayb->At);
        
    real_t wesum = 0.;
    int16_t * cl_intensities;
    NUC * cl_bases;
    const int ncpu = omp_get_max_threads();
    int th_id;
    int_fast32_t cl=0;
    MAT p[ncpu],V[ncpu];
    for( int i=0 ; i<ncpu ; i++){ p[i] = NULL; V[i] = NULL; }
    #pragma omp parallel for \
      default(shared) private(cl,cl_bases,cl_intensities,th_id) \
      reduction(+:wesum)
    for ( cl=0 ; cl<ncluster; cl++){
	th_id = omp_get_thread_num();
        cl_intensities = ayb->intensities.elt+cl*ncycle*NBASE;
        cl_bases = ayb->bases.elt + cl*ncycle; 
	p[th_id] =  processNew( AtLU, ayb->N, cl_intensities,p[th_id]);
        V[th_id] = accumulate_covariance(ayb->we->x[cl],p[th_id],ayb->lambda->x[cl],cl_bases,V[th_id]);
        wesum += ayb->we->x[cl];
    }
   
    // Accumulate per-thread results
    int lda = V[0]->nrow;
    MAT Vsum = new_MAT(lda,lda);
    for ( int i=0 ; i<ncpu ; i++){
	    if(NULL!=V[i]){
		    for ( int j=0 ; j<lda*lda ; j++){
			    Vsum->x[j] += V[i]->x[j];
		    }
	    }
    }
    for ( int i=0 ; i<ncpu ; i++){
    	free_MAT(p[i]);
    	free_MAT(V[i]);
    }
    free_MAT(AtLU.mat);
    free(AtLU.piv);

    // V is lower triangular
    for ( int i=0 ; i<lda ; i++){
       for ( int j=0 ; j<i ; j++){
	  Vsum->x[i*lda+j] = Vsum->x[j*lda+i];
       }
    }

    // Scale sum of squares to make covariance
    for ( uint32_t i=0 ; i<lda*lda ; i++){
            Vsum->x[i] /= wesum;
        }

    return Vsum;
}


