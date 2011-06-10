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

#include <math.h>
#include <omp.h>
#include "ayb.h"
#include "call_bases.h"
#include "estimate_lambda.h"
#include "process_intensities.h"
#include "estimateMPC.h"
#include "statistics.h"
#include "conjugate.h"


#define AYB_NITER   20


/*real_t initial_crosstalk[] = {
    2.0114300, 1.7217841, 0.06436576, 0.1126401,
    0.6919319, 1.8022413, 0.06436576, 0.0804572,
    0.2735545, 0.2252802, 1.39995531, 0.9976693,
    0.2896459, 0.2413716, 0.11264008, 1.3194981
};*/
// Bpert s_2_0050
/*real_t initial_crosstalk[] = {
   1.60, 1.23, -0.02, -0.03,
   0.27,  1.29, -0.03, -0.05,
   -0.10, -0.14, 0.86, 0.56,
   -0.10, -0.14, -0.02, 0.67
};*/
// HiSeq
real_t initial_crosstalk[] = {
     1.08, 1.13, 0.02, 0.01,
     0.20, 0.93, 0.02, 0.02,
     0.01, 0.02, 1.00, 0.53,
     0.00, 0.01, 0.05, 1.32
};

MAT initial_M = NULL;
MAT initial_At = NULL;

/* Basic functions */
AYB new_AYB(const uint32_t ncycle, const uint32_t ncluster){
    AYB ayb = malloc(sizeof(*ayb));
    if(NULL==ayb){ return NULL;}
    ayb->ncycle = ncycle;
    ayb->ncluster = ncluster;
    ayb->intensities = new_ARRAY(int16_t)(ncluster*ncycle*NBASE);
    ayb->bases = new_ARRAY(NUC)(ncluster*ncycle);
    ayb->quals = new_ARRAY(PHREDCHAR)(ncluster*ncycle);
    ayb->M = new_MAT(NBASE,NBASE);
    ayb->P = new_MAT(ncycle,ncycle);
    ayb->N = new_MAT(NBASE,ncycle);
    ayb->At = new_MAT(NBASE*ncycle,NBASE*ncycle);
    ayb->lambda = new_MAT(ncluster,1);
    ayb->we = new_MAT(ncluster,1);
    ayb->cycle_var = new_MAT(4*ncycle,1);

    // Filtering
    ayb->filtered = false;
    ayb->passed_filter = calloc(ncluster,sizeof(bool));

    ayb->index = 0;
    ayb->readnum = 0;    
    
    if( NULL==ayb->intensities.elt || NULL==ayb->bases.elt || NULL==ayb->M ||
        NULL==ayb->P || NULL==ayb->N || NULL==ayb->lambda  || NULL==ayb->passed_filter || NULL==ayb->At){
        goto cleanup;
    }
    
    ayb->coordinates = NULL;

    return ayb;
    
cleanup:
    free_AYB(ayb);
    return NULL;
}

void free_AYB(AYB ayb){
    if(NULL==ayb){ return;}
    free_ARRAY(int16_t)(ayb->intensities);
    free_ARRAY(NUC)(ayb->bases);
    free_ARRAY(PHREDCHAR)(ayb->quals);
    free_MAT(ayb->M);
    free_MAT(ayb->P);
    free_MAT(ayb->N);
    free_MAT(ayb->At);
    free_MAT(ayb->we);
    free_MAT(ayb->cycle_var);
    free_MAT(ayb->lambda);
    xfree(ayb->passed_filter);
    free_COORD(ayb->coordinates);
    xfree(ayb);
}

AYB copy_AYB(const AYB ayb){
    if(NULL==ayb){return NULL;}
    AYB ayb_copy = malloc(sizeof(*ayb));
    if(NULL==ayb_copy){ return NULL;}

    ayb_copy->ncycle = ayb->ncycle;
    ayb_copy->ncluster = ayb->ncluster;
    
    ayb_copy->intensities = copy_ARRAY(int16_t)(ayb->intensities);
    if(NULL!=ayb->intensities.elt && NULL==ayb_copy->intensities.elt){ goto cleanup;}

    ayb_copy->bases = copy_ARRAY(NUC)(ayb->bases);
    if(NULL!=ayb->bases.elt && NULL==ayb_copy->bases.elt){ goto cleanup;}
    
    ayb_copy->quals = copy_ARRAY(PHREDCHAR)(ayb->quals);
    if(NULL!=ayb->quals.elt && NULL==ayb_copy->quals.elt){ goto cleanup;}

    ayb_copy->M = copy_MAT(ayb->M);
    if(NULL!=ayb->M && NULL==ayb_copy->M){ goto cleanup;}
    
    ayb_copy->P = copy_MAT(ayb->P);
    if(NULL!=ayb->P && NULL==ayb_copy->P){ goto cleanup;}

    ayb_copy->N = copy_MAT(ayb->N);
    if(NULL!=ayb->N && NULL==ayb_copy->N){ goto cleanup;}
    
    ayb_copy->At = copy_MAT(ayb->At);
    if(NULL!=ayb->At && NULL==ayb_copy->At){ goto cleanup;}
    
    ayb_copy->we = copy_MAT(ayb->we);
    if(NULL!=ayb->we && NULL==ayb_copy->we){ goto cleanup;}
    
    ayb_copy->cycle_var = copy_MAT(ayb->cycle_var);
    if(NULL!=ayb->cycle_var && NULL==ayb_copy->cycle_var){ goto cleanup;}
    
    ayb_copy->lambda = copy_MAT(ayb->lambda);
    if(NULL!=ayb->lambda && NULL==ayb_copy->lambda){ goto cleanup;}
    
    ayb_copy->filtered = ayb->filtered;
    ayb_copy->passed_filter = calloc(ayb->ncluster,sizeof(bool));
    if(NULL==ayb_copy->passed_filter){ goto cleanup;}
    for ( int32_t i=0 ; i<ayb->ncluster ; i++){
        ayb_copy->passed_filter[i] = ayb->passed_filter[i];
    }
    
    ayb_copy->coordinates = copy_COORD(ayb->coordinates);
    if(NULL==ayb_copy->coordinates && NULL==ayb->coordinates){ goto cleanup;}
    
    ayb_copy->index = ayb->index;
    ayb_copy->readnum = ayb->readnum;
    
    return ayb_copy;

cleanup:
    free_AYB(ayb_copy);
    return NULL;
}

void show_AYB(XFILE * fp, const AYB ayb){
    validate(NULL!=fp,);
    validate(NULL!=ayb,);
    xfprintf(xstderr,"%u cycles from %u clusters\n",ayb->ncycle,ayb->ncluster);
    xfputs("M:\n",xstderr); show_MAT(xstderr,ayb->M,NBASE,NBASE);
    xfputs("P:\n",xstderr); show_MAT(xstderr,ayb->P,ayb->ncycle,ayb->ncycle);
    xfputs("N:\n",xstderr); show_MAT(xstderr,ayb->N,NBASE,8);
    xfputs("we:\n",xstderr); show_MAT(xstderr,ayb->we,8,1);
    xfputs("cycle_var:\n",xstderr); show_MAT(xstderr,ayb->cycle_var,8,1);
    xfputs("lambda:\n",xstderr); show_MAT(xstderr,ayb->lambda,8,1);
    xfputs("Bases:\n",xstderr); show_ARRAY(NUC)(stderr,ayb->bases,"",ayb->ncycle);
    xfputc('\n',xstderr);
    xfputs("Quality:\n",xstderr); show_ARRAY(PHREDCHAR)(stderr,ayb->quals,"",ayb->ncycle);
    xfputc('\n',xstderr);
    xfputs("Intensities:\n",xstderr); show_ARRAY(int16_t)(stderr,ayb->intensities," ",NBASE*ayb->ncycle);
    xfputc('\n',xstderr);
}
    
/* More complex initialisation */
AYB initialise_AYB(const CIFDATA cif){
    AYB ayb = new_AYB(cif_get_ncycle(cif),cif_get_ncluster(cif));
    if(NULL==ayb){ return NULL;}
    if(NULL==initial_M){ initial_M = new_MAT_from_array(NBASE,NBASE,initial_crosstalk);}
    
    timestamp("Transposing intensities\n",stderr);
    /*  Intensities in CIF are transpose of AYB format */
    encInt cifint = cif_get_const_intensities(cif);
    for ( uint32_t cy=0 ; cy<ayb->ncycle ; cy++){
        for ( uint32_t base=0 ; base<NBASE ; base++){
            for ( uint32_t cl=0 ; cl<ayb->ncluster ; cl++){
                ayb->intensities.elt[cl*ayb->ncycle*NBASE+cy*NBASE+base] = cifint.i16[cy*NBASE*ayb->ncluster+base*ayb->ncluster+cl];                
            }
        }
    }
    
    timestamp("Setting matrices\n",stderr);
    /* Initial cross-talk from default array */
    copyinto_MAT(ayb->M,initial_M);
    initial_At = new_MAT(4*ayb->ncycle,4*ayb->ncycle);
    const int lda = ayb->ncycle*4;
    for ( int i=0 ; i<ayb->ncycle ; i++){
	    for ( int j=0 ; j<4 ; j++){
		    for ( int k=0 ; k<4 ; k++){
	    		  initial_At->x[(i*4+j)*lda + i*4+k] = initial_M->x[j*4+k];
		    }
	    }
    }
    transpose_inplace(initial_At);
    copyinto_MAT(ayb->At,initial_At);
    /* Initial noise is zero */
    set_MAT(ayb->N,0.);
    /* Initial phasing */
    #warning "Phasing not properly initialised yet"
    ayb->P = identity_MAT(ayb->ncycle);
    /* Initial weights and cycles are all equal. Arbitrarily one */
    set_MAT(ayb->we,1.);
    set_MAT(ayb->cycle_var,1.0);
    
    /*  Call bases and lambda for each cluster */
    //MAT invAt = invert(ayb->At); //transpose(ayb->A);
    struct structLU AtLU = LUdecomposition(ayb->At);
    timestamp("Processing clusters\n",stderr);
    int_fast32_t cl,cy;
    NUC * cl_bases = NULL;
    PHREDCHAR * cl_quals = NULL;
    int th_id;
    int ncpu = omp_get_max_threads();
    MAT pcl_int[ncpu];
    for ( int i=0 ; i<ncpu ; i++){ pcl_int[i] = NULL; }
    #pragma omp parallel for \
      default(shared) private(cl,th_id,cl_bases,cl_quals,cy)
    for ( cl=0 ; cl<ayb->ncluster ; cl++){
	th_id = omp_get_thread_num();
	pcl_int[th_id] =  processNew( AtLU, ayb->N, ayb->intensities.elt+cl*ayb->ncycle*NBASE, pcl_int[th_id]);
        cl_bases = ayb->bases.elt + cl*ayb->ncycle;
        cl_quals = ayb->quals.elt + cl*ayb->ncycle;
        for ( cy=0 ; cy<ayb->ncycle ; cy++){
            cl_bases[cy] = call_base_simple(pcl_int[th_id]->x+cy*NBASE);
            cl_quals[cy] = 33;
        }
        //ayb->lambda->x[cl] = estimate_lambdaOLS(pcl_int,cl_bases);
	ayb->lambda->x[cl] = estimate_lambda_A ( ayb->intensities.elt+cl*ayb->ncycle*NBASE, ayb->N,ayb->At, cl_bases,  ayb->ncycle);
    }
    free_MAT(AtLU.mat);
    free(AtLU.piv);
    for ( int i=0 ; i<ncpu ; i++){
    	free_MAT(pcl_int[i]);
    }
    
    return ayb;
}

real_t update_cluster_weights(AYB ayb){
    validate(NULL!=ayb,NAN);
    const uint32_t ncluster = ayb->ncluster;
    const uint32_t ncycle   = ayb->ncycle;
    real_t sumLSS = 0.;
    
    //MAT invAt = invert(ayb->At); //transpose(ayb->A);
    struct structLU AtLU = LUdecomposition(ayb->At);
    int_fast32_t cl,i,idx;
    NUC * cycle_bases;
    int16_t * cycle_ints;
    const int ncpu = omp_get_max_threads();
    int th_id;
    MAT e[ncpu];
    for ( int i=0 ; i<ncpu ; i++){ e[i] = NULL; }
    /*  Calculate least squares error, using ayb->we as temporary storage */
    #pragma omp parallel for \
      default(shared) private(th_id,cl,cycle_ints,cycle_bases,i,idx) \
      reduction(+:sumLSS)
    for ( cl=0 ; cl<ncluster ; cl++){
	th_id = omp_get_thread_num();
        ayb->we->x[cl] = 0.;
        cycle_ints  = ayb->intensities.elt + cl*ncycle*NBASE;
        cycle_bases = ayb->bases.elt + cl*ncycle;
        e[th_id] = processNew( AtLU, ayb->N, cycle_ints, e[th_id]);
        for ( i=0 ; i<ncycle ; i++){
                e[th_id]->x[i*NBASE+cycle_bases[i]] -= ayb->lambda->x[cl];
        }
        for( idx=0 ; idx<NBASE*ncycle ; idx++){
            ayb->we->x[cl] += e[th_id]->x[idx]*e[th_id]->x[idx];
        }
        sumLSS += ayb->we->x[cl];
    }
    free_MAT(AtLU.mat);
    free(AtLU.piv);
    for ( int i=0 ; i<ncpu ; i++){
    	free_MAT(e[i]);
    }

    /* Calculate weight for each cluster */
    real_t meanLSSi = mean(ayb->we->x,ncluster);
    real_t varLSSi = variance(ayb->we->x,ncluster);
    for ( uint_fast32_t cl=0 ; cl<ncluster ; cl++){
        const real_t d = ayb->we->x[cl]-meanLSSi;
        ayb->we->x[cl] = cauchy(d*d,varLSSi);
    }
    //xfputs("Cluster weights:\n",xstderr);
    //show_MAT(xstderr,ayb->we,8,1);
    return sumLSS;
}

/*  Code for parameter estimation loop
 *  Equivalent to parameter_estimation_loop in old AYB
 *  Updates: M, P, N.
 *  Recalculates weights
 *  Scales all lambda by factor
 */
const real_t ridgeConst = 1000000.0;
MAT gblOmega = NULL;
real_t estimate_MPC( AYB ayb ){
    validate(NULL!=ayb,NAN);
    const uint32_t ncycle = ayb->ncycle;
    /*  Calculate new weights */
    //timestamp("Updating weights\n",stderr);
    real_t sumLSS = update_cluster_weights(ayb);
    
    /*  Precalculate terms for iteration */
    //timestamp("Calculating matrices\n",stderr);
    //timestamp("J\t",stderr);
    MAT J = calculateNewJ(ayb->lambda,ayb->bases,ayb->we,ncycle,NULL);
    //timestamp("K\t",stderr);
    MAT K = calculateNewK(ayb->lambda,ayb->bases,ayb->intensities,ayb->we,ncycle,NULL);
    //timestamp("Others\n",stderr);
    MAT Sbar = calculateSbar(ayb->lambda,ayb->we,ayb->bases,ncycle,NULL);
    MAT Ibar = calculateIbar(ayb->intensities,ayb->we,NULL);
    real_t Wbar = calculateWbar(ayb->we);
    real_t lambdaf = 1.;
    real_t * tmp = calloc(ncycle*ncycle*NBASE*NBASE,sizeof(real_t));
  
    MAT lhs = calculateLhs(Wbar, J, Sbar, NULL); 
    MAT rhs = calculateRhs( K, Ibar, NULL);
    for ( int i=0 ; i<lhs->nrow ; i++){
	    lhs->x[i*lhs->nrow+i] += ridgeConst;
    }
    for ( int i=0 ; i<initial_At->ncol ; i++){
        for ( int j=0 ; j<initial_At->nrow ; j++){
	    rhs->x[i*rhs->nrow+j] += ridgeConst * initial_At->x[i*initial_At->nrow+j];
	}
    }

    solver(lhs,rhs);
    for( int i=0 ; i<rhs->ncol ; i++){
        for( int j=0 ;j<rhs->ncol ; j++){
	    ayb->At->x[i*ayb->At->nrow+j] = rhs->x[i*rhs->nrow+j];
	}
	ayb->N->x[i] = rhs->x[i*rhs->nrow+rhs->nrow-1];
    }

    //real_t delta = calculateDeltaLSE( matMt, matP, matN, J, K, tmp);

    real_t det = normalise_MAT(ayb->At,3e-8);
    lambdaf *= det;
    // Scale lambdas by factor 
    scale_MAT(ayb->lambda,lambdaf);
    
    // Clean-up memory
    free_MAT(rhs);
    free_MAT(lhs);
    xfree(tmp);
    free_MAT(Ibar);
    free_MAT(Sbar);
    free_MAT(K);
    free_MAT(J);
    
    //xfprintf(xstderr,"Initial %e\tImprovement %e\t = %e\n",sumLSS,delta,sumLSS-delta);
    //xfprintf(xstderr,"Updated weights %e\n", update_cluster_weights(ayb));
    return sumLSS;//-delta;
}

real_t estimate_Bases(AYB ayb){
    validate(NULL!=ayb,NAN);
    //initialise_calibration();
    const uint32_t ncycle   = ayb->ncycle;
    const uint32_t ncluster = ayb->ncluster;
//timestamp("Calculating covariance\n",stderr);
    MAT V = calculate_covariance(ayb);
    // Scale is variance of residuals. Get from V matrices.
    for ( uint32_t i=0 ; i<4*ncycle ; i++){
        ayb->cycle_var->x[i] = V->x[i*ncycle*NBASE+i];
    }

    gblOmega = fit_omega(V,gblOmega,false);
    free_MAT(V); V = NULL;

    instrument( XFILE * matfp = xfopen("information.txt",XFILE_RAW,"w");
	        show_MAT(matfp,gblOmega,0,0);
                xfclose(matfp);
	  );

    //MAT invAt = invert(ayb->At); //transpose(ayb->A);
    struct structLU AtLU = LUdecomposition(ayb->At);

//timestamp("Base calling loop\n",stderr);
    //real_t qual[ncycle];
    int th_id;
    const int ncpu = omp_get_max_threads();
    int_fast32_t cl=0;
    NUC * bases = NULL;
    PHREDCHAR * phred = NULL;
    MAT pcl_int[ncpu];
    for ( int i=0 ; i<ncpu ; i++){ pcl_int[i] = NULL; }
    FILE * fp = fopen("processed.txt","w");
    #pragma omp parallel for \
      default(shared) private(cl,bases,phred,th_id)
    for ( cl=0 ; cl<ncluster ; cl++){
	th_id = omp_get_thread_num();
        bases = ayb->bases.elt + cl*ncycle;
        phred = ayb->quals.elt + cl*ncycle;
	pcl_int[th_id] =  processNew( AtLU, ayb->N, ayb->intensities.elt+cl*ayb->ncycle*NBASE, pcl_int[th_id]);
	ayb->lambda->x[cl] = estimate_lambda_A ( ayb->intensities.elt+cl*ayb->ncycle*NBASE, ayb->N, ayb->At, bases,  ayb->ncycle);
	// Call bases
	call_base(pcl_int[th_id],ayb->lambda->x[cl],gblOmega,bases);
	//for(int i=0 ; i<ncycle ; i++){ qual[i] = 40.0;}

	ayb->lambda->x[cl] = estimate_lambda_A ( ayb->intensities.elt+cl*ayb->ncycle*NBASE, ayb->N, ayb->At, bases,  ayb->ncycle);
	fprintf(fp,"%e %e %e %e\n",pcl_int[th_id]->x[0],pcl_int[th_id]->x[1],pcl_int[th_id]->x[2],pcl_int[th_id]->x[3]);
    }   
    fclose(fp);
//timestamp("Finished base calling\n",stderr);
    
    for ( int i=0 ; i<ncpu ; i++){ free_MAT(pcl_int[i]); }
    free_MAT(AtLU.mat);
    free(AtLU.piv);
    
    return NAN;
}
