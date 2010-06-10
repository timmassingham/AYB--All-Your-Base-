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

/*
 * Call base from processed intensities using minimum Least Squares.
 *  p       Processed intensities for given cycle
 *  lambda  Cluster brightness
 *  omega   Cycle specific inverse covariance matrix
 */
struct basequal call_base( const real_t * restrict p, const real_t lambda, const MAT omega){
    assert(NULL!=p);
    assert(NULL!=omega);
    assert(NBASE==omega->nrow);
    assert(NBASE==omega->ncol);
    extern AYBOPT aybopt;
    
    if(0==lambda){
        // Return random base with low quality
        int base = (int)(random()%4);
        struct basequal b = {base,33};
        return b;
    }
    
    int call = 0;
    real_t minstat = HUGE_VAL;
    real_t stat[NBASE] = { 0.,0.,0.,0.};
    for ( int i=0 ; i<NBASE ; i++){
        stat[i] = lambda * omega->x[i*NBASE+i];
        for ( int j=0 ; j<NBASE ; j++){
            stat[i] -= 2.0 * p[j] * omega->x[i*NBASE+j];
        }

        if(stat[i]<minstat){
            minstat = stat[i];
            call = i;
        }
    }
    /* Summation of probabilities for normalisation, 
     * having removed factor exp(-0.5*(K+lambda*minstat))
     */
    real_t tot = 0.;
    for ( int i=0 ; i<NBASE ; i++){
        tot += exp(-0.5*lambda*(stat[i]-minstat));
    }

    real_t K = xMy(p,omega,p);
    real_t maxprob = exp(-0.5*(K+lambda*minstat));

    /* Calculate posterior probability in numerically stable fashion
     * Note that maxp can be extremely small.
     */
    real_t post_prob = (maxprob<aybopt.mu) ?
       // Case probabilities small compared to mu
       (aybopt.mu + maxprob ) / (4.0*aybopt.mu + maxprob*tot) :
       // Case probabilities large compared to mu
       (aybopt.mu/maxprob + 1.) / (4.0*aybopt.mu/maxprob + tot);
    
    struct basequal b = {call,phredchar_from_prob(post_prob)};
    return b;
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
MAT * accumulate_covariance( const real_t we, const MAT p, const real_t lambda, const NUC * base, MAT * V){
    validate(NULL!=p,NULL);
    validate(NBASE==p->nrow,NULL);
    validate(lambda>=0.,NULL);
    const uint32_t ncycle = p->ncol;
    // Allocate memory for V if necessary
    if ( NULL==V ){
        V = calloc(ncycle+1,sizeof(*V));
        if(NULL==V){ return NULL;}
        for ( uint32_t cycle=0 ; cycle<ncycle ; cycle++){
            V[cycle] = new_MAT(NBASE,NBASE);
            if(NULL==V[cycle]){ goto cleanup;}
        }
    }
    // Perform accululation. V += R R^t
    // Note: R = P - \lambda I_b, where I_b is unit vector with b'th elt = 1
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        p->x[cy*NBASE+base[cy]] -= lambda;
    }
    for ( uint32_t cycle=0 ; cycle<ncycle ; cycle++){
        const int cybase = base[cycle];
        validate(cybase<NBASE,NULL);
        // P P^t
        for ( uint32_t i=0 ; i<NBASE ; i++){
            for ( uint32_t j=0 ; j<NBASE ; j++){
                V[cycle]->x[i*NBASE+j] += we * p->x[cycle*NBASE+i] * p->x[cycle*NBASE+j];
            }
        }
    }

    return V;
    
cleanup:
    for ( uint32_t cycle=0 ; cycle<ncycle ; cycle++){
        free_MAT(V[cycle]);
    }
    xfree(V);
    return NULL;
}


/*
 *  Calculate covariance of (processed residuals) 
 */
MAT * calculate_covariance( const AYB ayb){
    const uint32_t ncluster = ayb->ncluster;
    const uint32_t ncycle = ayb->ncycle;
    
    MAT * V = NULL;
    MAT p = NULL;
    // Create t(inverse(M)) and t(inverse(P))
    MAT Minv_t = transpose_inplace(invert(ayb->M));
    MAT Pinv_t = transpose_inplace(invert(ayb->P));
        
    real_t wesum = 0.;
    for ( uint32_t cl=0 ; cl<ncluster; cl++){
        const int16_t * cl_intensities = ayb->intensities.elt+cl*ncycle*NBASE;
        const NUC * cl_bases = ayb->bases.elt + cl*ncycle; 
        p = process_intensities(cl_intensities,Minv_t,Pinv_t,ayb->N,p);
        validate(NULL!=p,NULL);
        V = accumulate_covariance(ayb->we->x[cl],p,ayb->lambda->x[cl],cl_bases,V);
        validate(NULL!=V,NULL);
        wesum += ayb->we->x[cl];
    }
    free_MAT(p);
    free_MAT(Pinv_t);
    free_MAT(Minv_t);

    // Scale sum of squares to make covariance
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for ( uint32_t i=0 ; i<NBASE*NBASE ; i++){
            V[cy]->x[i] /= wesum;
        }
    }
    return V;
}


