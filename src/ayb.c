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
#include "ayb.h"
#include "call_bases.h"
#include "estimate_lambda.h"
#include "process_intensities.h"
#include "estimateMPC.h"
#include "statistics.h"


#define AYB_NITER   20


real_t initial_crosstalk[] = {
    2.0114300, 1.7217841, 0.06436576, 0.1126401,
    0.6919319, 1.8022413, 0.06436576, 0.0804572,
    0.2735545, 0.2252802, 1.39995531, 0.9976693,
    0.2896459, 0.2413716, 0.11264008, 1.3194981
};
MAT initial_M = NULL;


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
    ayb->lambda = new_MAT(ncluster,1);
    ayb->we = new_MAT(ncluster,1);
    ayb->cycle_var = new_MAT(ncycle,1);

    // Filtering
    ayb->filtered = false;
    ayb->passed_filter = calloc(ncluster,sizeof(bool));

    ayb->index = 0;
    ayb->readnum = 0;    
    
    if( NULL==ayb->intensities.elt || NULL==ayb->bases.elt || NULL==ayb->M ||
        NULL==ayb->P || NULL==ayb->N || NULL==ayb->lambda  || NULL==ayb->passed_filter){
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
    /* Initial noise is zero */
    set_MAT(ayb->N,0.);
    /* Initial phasing */
    #warning "Phasing not properly initialised yet"
    ayb->P = identity_MAT(ayb->ncycle);
    /* Initial weights and cycles are all equal. Arbitrarily one */
    set_MAT(ayb->we,1.);
    set_MAT(ayb->cycle_var,1);
    
    /*  Call bases and lambda for each cluster */
    MAT pcl_int = NULL; // Shell for processed intensities
    MAT Minv_t = transpose_inplace(invert(ayb->M));
    MAT Pinv_t = transpose_inplace(invert(ayb->P));
    timestamp("Processing clusters\n",stderr);
    for ( uint32_t cl=0 ; cl<ayb->ncluster ; cl++){
        pcl_int = process_intensities(ayb->intensities.elt+cl*ayb->ncycle*NBASE,Minv_t,Pinv_t,ayb->N,pcl_int);
        NUC * cl_bases = ayb->bases.elt + cl*ayb->ncycle;
        PHREDCHAR * cl_quals = ayb->quals.elt + cl*ayb->ncycle;
        for ( uint32_t cy=0 ; cy<ayb->ncycle ; cy++){
            cl_bases[cy] = call_base_simple(pcl_int->x+cy*NBASE);
            cl_quals[cy] = 33;
        }
        ayb->lambda->x[cl] = estimate_lambdaOLS(pcl_int,cl_bases);
    }
    free_MAT(pcl_int);
    free_MAT(Pinv_t);
    free_MAT(Minv_t);
    
    return ayb;
}

real_t update_cluster_weights(AYB ayb){
    validate(NULL!=ayb,NAN);
    const uint32_t ncluster = ayb->ncluster;
    const uint32_t ncycle   = ayb->ncycle;
    real_t sumLSS = 0.;
    
    MAT e = NULL;
    /*  Calculate least squares error, using ayb->we as temporary storage */
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        ayb->we->x[cl] = 0.;
        int16_t * cycle_ints  = ayb->intensities.elt + cl*ncycle*NBASE;
        NUC *     cycle_bases = ayb->bases.elt + cl*ncycle;
        e = expected_intensities(ayb->lambda->x[cl],cycle_bases,ayb->M,ayb->P,ayb->N,e);
        for( uint32_t idx=0 ; idx<NBASE*ncycle ; idx++){
            real_t tmp = cycle_ints[idx] - e->x[idx];
            ayb->we->x[cl] += tmp*tmp;
        }
        sumLSS += ayb->we->x[cl];
    }
    free_MAT(e);
    
    /* Calculate weight for each cluster */
    real_t meanLSSi = mean(ayb->we->x,ncluster);
    real_t varLSSi = variance(ayb->we->x,ncluster);
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
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
real_t estimate_MPC( AYB ayb ){
    validate(NULL!=ayb,NAN);
    const uint32_t ncycle = ayb->ncycle;
    /*  Calculate new weights */
    //timestamp("Updating weights\n",stderr);
    real_t sumLSS = update_cluster_weights(ayb);
    
    /*  Precalculate terms for iteration */
    //timestamp("Calculating matrices\n",stderr);
    //timestamp("J\t",stderr);
    MAT J = calculateJ(ayb->lambda,ayb->we,ayb->bases,ncycle,NULL);
    MAT Jt = transpose(J);
    //timestamp("K\t",stderr);
    MAT K = calculateK(ayb->lambda,ayb->we,ayb->bases,ayb->intensities,ncycle,NULL);
    //timestamp("Others\n",stderr);
    MAT Kt = transpose(K);
    MAT Sbar = calculateSbar(ayb->lambda,ayb->we,ayb->bases,ncycle,NULL);
    MAT SbarT = transpose(Sbar);
    MAT Ibar = calculateIbar(ayb->intensities,ayb->we,NULL);
    MAT IbarT = transpose(Ibar);
    real_t Wbar = calculateWbar(ayb->we);
    real_t lambdaf = 1.;
    real_t * tmp = calloc(ncycle*ncycle*NBASE*NBASE,sizeof(real_t));
    /* Convenience terms for M, P and N */
    MAT matMt = transpose_inplace(ayb->M);
    MAT matP = ayb->P;
    MAT matN = ayb->N;
    
    //xfputs("Staring main loop",xstderr);
    /*  Main iteration */
    MAT mlhs=NULL,mrhs=NULL,plhs=NULL,prhs=NULL;
    const uint32_t lda = NBASE + ncycle;
    real_t det=0.;
//xfprintf(xstderr,"Starting\tdelta = %e\n",calculateDeltaLSE( matMt, matP, matN, J, K, tmp));
    for( uint32_t i=0 ; i<AYB_NITER ; i++){
        /*
         *  Solution for phasing and constant noise
         */
        plhs = calculatePlhs(Wbar,Sbar,matMt,J,tmp,plhs);
        prhs = calculatePrhs(Ibar,matMt,K,tmp,prhs);
        solverSVD(plhs,prhs,tmp);
        for ( uint32_t i=0 ; i<ncycle ; i++){
            for(uint32_t j=0 ; j<ncycle ; j++){
                matP->x[i*ncycle+j] = (prhs->x[i*ncycle+j]>=0.)?prhs->x[i*ncycle+j]:0.;
            }
        }

        // Scaling so det(P) = 1
        det = normalise_MAT(matP,3e-8);
        scale_MAT(J,det*det); scale_MAT(Jt,det*det);
        scale_MAT(K,det);     scale_MAT(Kt,det);
        scale_MAT(Sbar,det);  scale_MAT(SbarT,det);
        lambdaf *= det;
        
        /*
         *  Solution for cross-talk and constant noise
         */
        mlhs = calculateMlhs(ayb->cycle_var,Wbar,SbarT,matP,Jt,tmp,mlhs);
        mrhs = calculateMrhs(ayb->cycle_var,IbarT,matP,Kt,tmp,mrhs);
        solverSVD(mlhs,mrhs,tmp);

        for(uint32_t i=0 ; i<NBASE ; i++){
            for(uint32_t j=0 ; j<NBASE ; j++){
                matMt->x[i*NBASE+j] = mrhs->x[i*lda+j];
            }
        }
        for(uint32_t i=0 ; i<NBASE ; i++){
            for(uint32_t j=0 ; j<ncycle ; j++){
                matN->x[j*NBASE+i] = mrhs->x[i*lda+j+NBASE];
            }
        }
        // Scaling so det(M) = 1
        det = normalise_MAT(matMt,3e-8);
        scale_MAT(J,det*det); scale_MAT(Jt,det*det);
        scale_MAT(K,det);     scale_MAT(Kt,det);
        scale_MAT(Sbar,det);  scale_MAT(SbarT,det);
        lambdaf *= det;
        //xfprintf(xstderr,"... done %u\tdelta = %e\n",i,calculateDeltaLSE( matMt, matP, matN, J, K, tmp));
    }
    real_t delta = calculateDeltaLSE( matMt, matP, matN, J, K, tmp);

    // Transpose Mt back to normal form 
    matMt = transpose_inplace(matMt);
    // Scale lambdas by factor 
    scale_MAT(ayb->lambda,lambdaf);
    
    // Clean-up memory
    free_MAT(prhs);
    free_MAT(plhs);
    free_MAT(mrhs);
    free_MAT(mlhs);
    xfree(tmp);
    free_MAT(IbarT);
    free_MAT(Ibar);
    free_MAT(SbarT);
    free_MAT(Sbar);
    free_MAT(Kt);
    free_MAT(K);
    free_MAT(Jt);
    free_MAT(J);
    
    //xfprintf(xstderr,"Initial %e\tImprovement %e\t = %e\n",sumLSS,delta,sumLSS-delta);
    //xfprintf(xstderr,"Updated weights %e\n", update_cluster_weights(ayb));
    return sumLSS-delta;
}

real_t estimate_Bases(AYB ayb){
    validate(NULL!=ayb,NAN);
    const uint32_t ncycle   = ayb->ncycle;
    const uint32_t ncluster = ayb->ncluster;
//timestamp("Calculating covariance\n",stderr);
    MAT * V = calculate_covariance(ayb);
    // Scale is variance of residuals. Get from V matrices.
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        ayb->cycle_var->x[cy] = 0.;
        for ( uint32_t b=0 ; b<NBASE ; b++){
            ayb->cycle_var->x[cy] += V[cy]->x[b*NBASE+b];
        }
        ayb->cycle_var->x[cy] = ayb->cycle_var->x[cy];
    }

    // Invert variance matrices
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        MAT a = V[cy];
        V[cy] = invert(a);
        free_MAT(a);
    }

    MAT pcl_int = NULL; // Shell for processed intensities
    MAT Minv_t = transpose_inplace(invert(ayb->M));
    MAT Pinv_t = transpose_inplace(invert(ayb->P));

//timestamp("Base calling loop\n",stderr);
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        NUC * bases = ayb->bases.elt + cl*ncycle;
        PHREDCHAR * quals = ayb->quals.elt + cl*ncycle;
        pcl_int = process_intensities(ayb->intensities.elt+cl*ncycle*NBASE,Minv_t,Pinv_t,ayb->N,pcl_int);
        //ayb->lambda->x[cl] = estimate_lambdaGWLS(pcl_int,bases,ayb->lambda->x[cl],ayb->cycle_var->x,V);
        ayb->lambda->x[cl] = estimate_lambdaWLS(pcl_int,bases,ayb->lambda->x[cl],ayb->cycle_var->x);
        for ( uint32_t cy=0 ; cy<ncycle ; cy++){
            struct basequal bq = call_base(pcl_int->x+cy*NBASE,ayb->lambda->x[cl],V[cy]);
            bases[cy] = bq.base;
            quals[cy] = bq.qual;
        }
        //ayb->lambda->x[cl] = estimate_lambdaGWLS(pcl_int,bases,ayb->lambda->x[cl],ayb->cycle_var->x,V);
        ayb->lambda->x[cl] = estimate_lambdaWLS(pcl_int,bases,ayb->lambda->x[cl],ayb->cycle_var->x);
    }    
//timestamp("Finished base calling\n",stderr);
    
    free_MAT(pcl_int);
    free_MAT(Pinv_t);
    free_MAT(Minv_t);
    for(uint32_t cy=0 ; cy<ncycle ; cy++){
        free_MAT(V[cy]);
    }
    free(V);
    
    return NAN;
}
