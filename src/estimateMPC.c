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
#include <strings.h>
#include "ayb.h"
#include "estimateMPC.h"
#include "utility.h"
#include "matrix.h"
#include "statistics.h"
#include "nuc.h"
#include "lapack.h"


static inline real_t sum_squares(const real_t * x, const uint32_t n){
    validate(NULL!=x,NAN);
    real_t res = 0.;
    for ( uint32_t i=0 ; i<n ; i++){
        res += x[i] * x[i];
    }
    return res;
}


// Crude method of obtaining weights
MAT calculateWe( const MAT lssi, MAT we){
    validate(NULL!=lssi,NULL);
    const uint32_t ncluster = lssi->nrow;
    if(NULL==we){
        we = new_MAT(ncluster,1);
        validate(NULL!=we,NULL);
    }
    bzero(we->x,ncluster*sizeof(real_t));
    
    real_t meanLSSi = mean(lssi->x,ncluster);
    real_t varLSSi = variance(lssi->x,ncluster);
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        const real_t d = lssi->x[cl]-meanLSSi;
        we->x[cl] = cauchy(d*d,varLSSi);
    }
    return we;
}


MAT calculateIbar( const ARRAY(int16_t) intmat, const MAT we, MAT Ibar){
    validate(NULL!=intmat.elt,NULL);
    validate(NULL!=we,NULL);
    const uint32_t ncluster = we->nrow;

    const uint32_t ncycle = intmat.nelt/(ncluster*NBASE);
    if(NULL==Ibar){
        Ibar = new_MAT(NBASE,ncycle);
        validate(NULL!=Ibar,NULL);
    }
    validate(Ibar->nrow==NBASE,NULL);
    validate(Ibar->ncol==ncycle,NULL);
    bzero(Ibar->x,Ibar->nrow*Ibar->ncol*sizeof(real_t));
    
    for( uint32_t cl=0 ; cl<ncluster ; cl++){
        for( uint32_t idx=0 ; idx<NBASE*ncycle ; idx++){
            Ibar->x[idx] += intmat.elt[cl*NBASE*ncycle+idx] * we->x[cl];
        }
    }
    
    return Ibar;
}

MAT calculateSbar( const MAT lambda, const MAT we, const ARRAY(NUC) bases, const uint32_t ncycle, MAT Sbar){
    validate(NULL!=lambda,NULL);
    validate(NULL!=we,NULL);
    validate(NULL!=bases.elt,NULL);
    const uint32_t ncluster = lambda->nrow;
    validate(ncluster==we->nrow,NULL);
    if(NULL==Sbar){
        Sbar = new_MAT(NBASE,ncycle);
        validate(NULL!=Sbar,NULL);
    }
    validate(Sbar->nrow==NBASE,NULL);
    validate(Sbar->ncol==ncycle,NULL);
    bzero(Sbar->x,Sbar->nrow*Sbar->ncol*sizeof(real_t));
    
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        for ( uint32_t cy=0 ; cy<ncycle ; cy++){
            int base = bases.elt[cl*ncycle+cy];
            Sbar->x[cy*NBASE+base] += we->x[cl] * lambda->x[cl];
        }
    }
    
    return Sbar;
}
real_t calculateWbar( const MAT we){
    validate(NULL!=we,NAN);
    const uint32_t ncluster = we->nrow;
    
    real_t wbar = 0.;
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        wbar += we->x[cl];
    }
    return wbar;
}
        

MAT calculateJ( const MAT lambda, const MAT we, const ARRAY(NUC) bases, const uint32_t ncycle, MAT J){
    validate(NULL!=lambda,NULL);
    validate(NULL!=we,NULL);
    validate(NULL!=bases.elt,NULL);
    const uint32_t ncluster = lambda->nrow;

    if(NULL==J){
        J = new_MAT(NBASE*NBASE,ncycle*ncycle);
        validate(NULL!=J,NULL);
    }
    bzero(J->x,J->nrow*J->ncol*sizeof(real_t));
    
    const uint32_t lda = NBASE*NBASE;
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        const real_t welam = we->x[cl] * lambda->x[cl] * lambda->x[cl];
        for ( uint32_t cy=0 ; cy<ncycle ; cy++){
            const int base = bases.elt[cl*ncycle+cy];
            const int offset = cy*ncycle*NBASE*NBASE + base*NBASE;
            for ( uint32_t cy2=0 ; cy2<ncycle ; cy2++){
                const int base2 = bases.elt[cl*ncycle+cy2];
                J->x[offset+cy2*lda+base2] += welam;
            }
        }
    }
    
    return J;
}

// tmp should be NBASE*NBASE+ncycle*ncycle
real_t calculateDeltaLSE(const MAT Mt, const MAT P, const MAT N, const MAT J, const MAT K, real_t * tmp){
    // BLAS definitions
    const real_t alpha = 1.0;
    const real_t  beta = 0.0;
    // Sanity checks
    validate(NULL!=Mt,NAN);
    validate(NULL!=P,NAN);
    validate(NULL!=N,NAN);
    validate(NULL!=J,NAN);
    validate(NULL!=K,NAN);
    validate(NULL!=tmp,NAN);
    const int ncycle = P->nrow;
    const int nbase = NBASE;

    real_t delta = -sum_squares(N->x,ncycle*NBASE); // tr(N^tN) = 1^t (N o N) 1 = sum_squares
    delta -= 2.0 * xMy(Mt->x,K,P->x);
    // PP^t
    gemm(LAPACK_NOTRANS,LAPACK_TRANS,&ncycle,&ncycle,&ncycle,&alpha,P->x,&ncycle,P->x,&ncycle,&beta,tmp+NBASE*NBASE,&ncycle);
    // MtM
    gemm(LAPACK_NOTRANS,LAPACK_TRANS,&nbase, &nbase, &nbase, &alpha,Mt->x,&nbase, Mt->x,&nbase, &beta,tmp,&nbase);
    delta += xMy(tmp,J,tmp+NBASE*NBASE);
    
    return delta;
}

MAT calculateK( const MAT lambda, const MAT we, const ARRAY(NUC) bases, const ARRAY(int16_t) ints, const uint32_t ncycle, MAT K){
    validate(NULL!=lambda,NULL);
    validate(NULL!=we,NULL);
    validate(NULL!=ints.elt,NULL);
    const uint32_t ncluster = lambda->nrow;

    if(NULL==K){
        K = new_MAT(NBASE*NBASE,ncycle*ncycle);
        validate(NULL!=K,NULL);
    }
    bzero(K->x,K->nrow*K->ncol*sizeof(real_t));
    
    const uint32_t lda = NBASE*NBASE;
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        const real_t welam = we->x[cl] * lambda->x[cl];
        for ( uint32_t cy=0 ; cy<ncycle ; cy++){
            const uint32_t ioffset = cl*NBASE*ncycle + cy*NBASE;
            for ( uint32_t cy2=0 ; cy2<ncycle ; cy2++){
                const int base = bases.elt[cl*ncycle+cy2];
                const uint32_t koffset = cy*lda*ncycle + cy2*lda + base;
                for ( uint32_t ch=0 ; ch<NBASE ; ch++){
                    K->x[ koffset + ch*NBASE] += welam * ints.elt[ioffset + ch];
                }
            }
        }
    }
    
    return K;
}


// tmp should be NBASE*NBASE+ncycle*ncycle
MAT calculateMlhs( const MAT var, const real_t wbar, const MAT SbarT, const MAT P, const MAT Jt, real_t * tmp, MAT lhs){
    validate(NULL!=var,NULL);
    validate(SbarT!=NULL,NULL);
    validate(P!=NULL,NULL);
    validate(Jt!=NULL,NULL);
    validate(tmp!=NULL,NULL);
    const int unit = 1;
    const real_t alpha = 1.0;
    const real_t  beta = 0.0;
    const int ncycle = P->nrow;
    
    if(NULL==lhs){
        lhs = new_MAT(NBASE+ncycle,NBASE+ncycle);
        validate(NULL!=lhs,NULL);
    }
    bzero(lhs->x,lhs->nrow*lhs->ncol*sizeof(real_t));

    const uint32_t lda = NBASE+ncycle;    
    // Reshape Jvec(P diag(var) Pt) via P diag(sqrt(v)) %*% ( P diag(sqrt(v)) )^t
    // Scale P
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for (uint32_t cy2=0 ; cy2<ncycle ; cy2++){
            P->x[cy*ncycle+cy2] *= sqrt(var->x[cy]);
        }
    }
    gemm(LAPACK_NOTRANS,LAPACK_TRANS,&ncycle,&ncycle,&ncycle,&alpha,P->x,&ncycle,P->x,&ncycle,&beta,tmp+NBASE*NBASE,&ncycle);
    gemv(LAPACK_TRANS,&Jt->nrow,&Jt->ncol,&alpha,Jt->x,&Jt->nrow,tmp+NBASE*NBASE,&unit,&beta,tmp,&unit);
    for ( uint32_t cy=0 ; cy<NBASE ; cy++){
        for ( uint32_t base=0 ; base<NBASE ; base++){
            lhs->x[cy*lda+base] = tmp[cy*NBASE+base];
        }
    }
    // Unscale P
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for (uint32_t cy2=0 ; cy2<ncycle ; cy2++){
            P->x[cy*ncycle+cy2] /= sqrt(var->x[cy]);
        }
    }
    
    // Sbar %*% P %*% diag(var)
    gemm(LAPACK_TRANS,LAPACK_NOTRANS,&SbarT->ncol,&P->ncol,&P->nrow,&alpha,SbarT->x,&SbarT->nrow,P->x,&P->nrow,&beta,lhs->x+NBASE*lda,&lhs->nrow);
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for ( uint32_t base=0 ; base<NBASE ; base++){
            lhs->x[NBASE*lda+cy*lda+base] *= var->x[cy];
        }
    }
    // Copy Sbar %*% P into new bit of array
{
    const uint32_t offset1 = lda * NBASE;
    const uint32_t offset2 = NBASE;
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for ( uint32_t base=0 ; base<NBASE ; base++){
            lhs->x[offset2+base*lda+cy] = lhs->x[offset1+cy*lda+base];
        }
    }
}
    // Id matrix
{
    const uint32_t offset = NBASE*lda+NBASE;
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        lhs->x[offset+cy*lda+cy] = wbar*var->x[cy];
    }
}
    
    return lhs;
}


// tmp should be NBASE*NBASE long
MAT calculateMrhs( const MAT var, const MAT IbarT, const MAT P, const MAT Kt, real_t * tmp, MAT rhs){
    validate(NULL!=var,NULL);
    validate(NULL!=IbarT,NULL);
    validate(NULL!=P,NULL);
    validate(NULL!=Kt,NULL);
    validate(NULL!=tmp,NULL);
    const real_t alpha = 1.0;
    const real_t  beta = 0.0;

    const uint32_t ncycle = P->nrow;
    const int lda = NBASE + ncycle;
    if(NULL==rhs){
        rhs = new_MAT(lda,NBASE);
        validate(NULL!=rhs,NULL);
    }
    bzero(rhs->x,rhs->nrow*rhs->ncol*sizeof(real_t));
    
    // Reshape KVec(Pdiag(v))
    // Scale P
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for (uint32_t cy2=0 ; cy2<ncycle ; cy2++){
            P->x[cy*ncycle+cy2] *= var->x[cy];
        }
    }
    gemv(LAPACK_TRANS,&Kt->nrow,&Kt->ncol,&alpha,Kt->x,&Kt->nrow,P->x,LAPACK_UNIT,&beta,tmp,LAPACK_UNIT);
    for ( uint32_t base1=0 ; base1<NBASE ; base1++){
        for ( uint32_t base2=0 ; base2<NBASE ; base2++){
            rhs->x[base1*lda+base2] = tmp[base1*NBASE+base2];
        }
    }
    // Unscale P
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for (uint32_t cy2=0 ; cy2<ncycle ; cy2++){
            P->x[cy*ncycle+cy2] /= var->x[cy];
        }
    }
    
    // Copy in diag(v)*IbarT
    for ( uint32_t base=0 ; base<NBASE ; base++){
        for ( uint32_t cy=0 ; cy<ncycle ; cy++){
            rhs->x[NBASE+base*lda+cy] = IbarT->x[base*ncycle+cy] * var->x[cy];
        }
    }
    return rhs;
}


// tmp should be NBASE*NBASE+ncycle*ncycle
MAT calculatePlhs( const real_t wbar, const MAT Sbar, const MAT Mt, const MAT J, real_t * tmp, MAT lhs){
    validate(Sbar!=NULL,NULL);
    validate(Mt!=NULL,NULL);
    validate(J!=NULL,NULL);
    validate(tmp!=NULL,NULL);
    const real_t alpha = 1.0;
    const real_t  beta = 0.0;
    const int ncycle = Sbar->ncol;
    const int nbase = NBASE;
    
    if(NULL==lhs){
        //lhs = new_MAT(NBASE+ncycle,NBASE+ncycle);
	lhs = new_MAT(ncycle,ncycle);
        validate(NULL!=lhs,NULL);
    }
    bzero(lhs->x,lhs->nrow*lhs->ncol*sizeof(real_t));

    const int lda = NBASE+ncycle;    
    // Reshape Jtvec(MtM)
    gemm(LAPACK_NOTRANS, LAPACK_TRANS,&nbase, &nbase, &nbase, &alpha,Mt->x,&nbase, Mt->x,  &nbase, &beta, tmp+ncycle*ncycle,&nbase);
    gemv(LAPACK_TRANS,&J->nrow,&J->ncol,&alpha,J->x,&J->nrow,tmp+ncycle*ncycle,LAPACK_UNIT,&beta,tmp,LAPACK_UNIT);
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for ( uint32_t cy2=0 ; cy2<ncycle ; cy2++){
            lhs->x[cy*ncycle+cy2] = tmp[cy*ncycle+cy2];
        }
    }
/*    // M %*% Sbar
    gemm(LAPACK_TRANS,LAPACK_NOTRANS,&nbase,&ncycle,&nbase,&alpha,Mt->x,&nbase,Sbar->x,&nbase,&beta,lhs->x+ncycle,&lda);
    // Copy M %*% Sbar  into new bit of array
{
    const uint32_t offset1 = lda * ncycle;
    const uint32_t offset2 = ncycle;
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for ( uint32_t base=0 ; base<NBASE ; base++){
            lhs->x[offset1+base*lda+cy] = lhs->x[offset2+cy*lda+base];
        }
    }
}
    // Id matrix
{
    const uint32_t offset = ncycle*lda+ncycle;
    for ( uint32_t base=0 ; base<NBASE ; base++){
        lhs->x[offset+base*lda+base] = wbar;
    }
}*/
    
    return lhs;
}


// tmp should be ncycle*ncycle long
MAT calculatePrhs( const MAT Ibar, const MAT Mt, const MAT K, real_t * tmp, MAT rhs){
    validate(NULL!=Ibar,NULL);
    validate(NULL!=Mt,NULL);
    validate(NULL!=K,NULL);
    validate(NULL!=tmp,NULL);
    const real_t alpha = 1.0;
    const real_t  beta = 0.0;

    const uint32_t ncycle = Ibar->ncol;
    //const int lda = NBASE + ncycle;
    if(NULL==rhs){
        //rhs = new_MAT(lda,ncycle);
	rhs = new_MAT(ncycle,ncycle);
        validate(NULL!=rhs,NULL);
    }
    bzero(rhs->x,rhs->nrow*rhs->ncol*sizeof(real_t));
    
    // Reshape KtVec(M)
    gemv(LAPACK_TRANS,&K->nrow,&K->ncol,&alpha,K->x,&K->nrow,Mt->x,LAPACK_UNIT,&beta,tmp,LAPACK_UNIT);
    for ( uint32_t cy1=0 ; cy1<ncycle ; cy1++){
        for ( uint32_t cy2=0 ; cy2<ncycle ; cy2++){
            rhs->x[cy1*ncycle+cy2] = tmp[cy1*ncycle+cy2];
        }
    }
/*
    // Copy in Ibar
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
        for ( uint32_t base=0 ; base<NBASE ; base++){
            rhs->x[ncycle+cy*lda+base] = Ibar->x[cy*NBASE+base];
        }
    }*/
    return rhs;
}

/* Assumes that LHS is positive-definite, which problems being considered are.
 * Both lhs and rhs are destructively updated
 */
int solver( MAT lhs, MAT rhs){
    validate(NULL!=lhs,-4);
    validate(NULL!=rhs,-6);
    validate(lhs->nrow==lhs->ncol,-2);
    validate(lhs->ncol==rhs->nrow,-3);

    int info = 0;
    posv(LAPACK_UPPER, &lhs->ncol, &rhs->ncol, lhs->x, &lhs->nrow, rhs->x, &rhs->nrow, &info);
    if(0!=info){ fprintf(stderr,"Error in solver, info = %d\n",info);}
    return info;
}


// Solve linear system using SVD
// tmp should be 6*N
int solverSVD(MAT lhs, MAT rhs, real_t * tmp){
    validate(NULL!=lhs,-4);
    validate(NULL!=rhs,-6);
    validate(NULL!=tmp,-11);
    
    const int N = lhs->nrow;
    int INFO=0,RANK=0,IWORK=5*N;
    double RCOND = 3e-8;
    gelss(&lhs->nrow,&lhs->ncol,&rhs->ncol,lhs->x,&lhs->nrow,
                      rhs->x,&rhs->nrow,tmp,&RCOND,&RANK,tmp+N,&IWORK,&INFO);
    return INFO;
}

void dnnls_(double *A, int * MDA, int * M, int* N, double * B, double * X, double * RNORM, double * W, double * ZZ, int * INDEX, int * MODE);

int solveNonneg(MAT lhs, MAT rhs, real_t *tmp){
    const int N = lhs->nrow;
    real_t RNORM;
    real_t W[N];
    real_t ZZ[N];
    int INDEX[N],MODE;

    double * lhs_tmp = malloc(lhs->nrow*lhs->ncol*sizeof(double));
    double * rhs_tmp = malloc(rhs->nrow*sizeof(double));
    for( int cy=0 ; cy<N ; cy++){
        memcpy(lhs_tmp,lhs->x,lhs->nrow*lhs->ncol*sizeof(double));
        memcpy(rhs_tmp,rhs->x+cy*N,rhs->nrow*sizeof(double));
        dnnls_(lhs_tmp,&N,&N,&N,rhs_tmp,tmp+cy*N,&RNORM,W,ZZ,INDEX,&MODE);
    }
    fprintf(stdout,"Solution result = %d\n",MODE);
    free(lhs_tmp);
    free(rhs_tmp);
    return MODE;
}

#ifdef TEST
#include <stdio.h>
static  int bases[] = {4, 3, 4, 4, 4, 2, 4, 2, 4, 1, 2, 4, 4, 1, 3};
static  real_t we[] = {0.33081400, 0.29916646, 0.79687810, 0.04376649, 0.04898477};
static  real_t lambda[] = {2.852168, 11.890647, 10.675880, 5.542096, 3.245923};
static  real_t ints_t[] = { 
    0.32822142, 0.3611550, 0.4056395, 0.6217474, 0.3701723, 0.54075999, 0.9004344, 0.8360398, 0.7729722, 0.9753776, 0.06152802, 0.5637200,
    0.03651698, 0.8598354, 0.4986751, 0.7010750, 0.3627614, 0.06679484, 0.9104111, 0.8990147, 0.4424424, 0.2002516, 0.11924224, 0.2969899,
    0.54164990, 0.9133029, 0.5054887, 0.5377883, 0.2405506, 0.58718250, 0.3922185, 0.1585344, 0.8682236, 0.3153632, 0.86569760, 0.8310858,
    0.61345680, 0.9237610, 0.9891983, 0.1324787, 0.3127209, 0.69892800, 0.2603790, 0.2256760, 0.8400769, 0.9111117, 0.27118720, 0.5472079,
    0.94282040, 0.4265982, 0.1394416, 0.7856237, 0.2671884, 0.11008570, 0.4843887, 0.1912447, 0.7967254, 0.5434273, 0.83477990, 0.7742892
};
static real_t SbarT[] = {
    0.4360868, 0.1086585, 0.3703858,
    0.5885691, 0.4436660, 0.5181457,
    0.4142125, 0.4619443, 0.4768580,
    0.3178741, 0.3513695, 0.1786516
};
static real_t P[] = {
    0.9445303038, 0.70435005, 0.74991005,
    0.0006557733, 0.05396262, 0.05993442,
    0.0018982966, 0.08522693, 0.54507703
};
static real_t Ibar[] = {
    0.9058495, 0.3065012, 0.9675010, 0.5600422,
    0.1745864, 0.3782057, 0.1538957, 0.5697690,
    0.5722346, 0.1426533, 0.5951879, 0.8001268
};
static real_t M[] = {
     0.8390111, 0.8550323, 0.6927066, 0.7199632,
     0.4642314, 0.8833080, 0.1549395, 0.1361410,
     0.1612448, 0.3195291, 0.5993373, 0.4147634,
     0.8614796, 0.3521542, 0.2954812, 0.7250016
};
/* Result

J matrix:
1:     1.34     0.00     0.00     0.00     0.52     0.00     0.00     0.00     0.00
2:     0.00     1.34     0.00     0.00     0.00     0.00     0.00     0.00     0.00
3:     0.00     0.00     0.00     0.00     0.00     0.52     0.00     0.00     0.00
4:     0.00     0.00     1.34     0.52     0.00     0.00     0.00     0.00     0.00
5:     0.00     0.00     0.00     1.34     0.00     0.00     0.00     0.00     0.00
6:     0.00     0.00     0.00     0.00    92.17     0.00     0.00     0.00    42.30
7:     0.00     0.00     0.00     0.00     0.00     0.00     0.00     0.00     0.00
8:     0.00     0.00     0.00    90.82     0.00    92.17    42.30    42.30     0.00
9:     0.00     0.00     0.00     0.00     0.00     0.00     0.00     0.52     0.00
10:     0.00     0.00     0.00     0.00     0.00     0.00     0.00     0.00     0.00
11:     0.00     0.00     0.00     0.00     2.69     0.00     0.00     0.00     0.52
12:     0.00     0.00     0.00     2.69     0.00     2.69     0.52     0.00     0.00
13:     0.00     0.52     0.00     0.00     0.00     0.00     1.34     0.00     0.00
14:     0.00    90.82    42.30     0.00     0.00    42.30     0.00    92.17     0.00
15:     0.00     2.69     0.52     0.00     0.00     0.00     0.00     2.69     0.00
16:   136.33    42.30    93.51    42.30    42.30     0.00    93.51     0.00    94.86

K matrix:
1:     0.15     0.15     0.00     0.08     0.04     0.00     0.20     0.13     0.00
2:     0.00     4.76     0.13     0.00     2.12     1.29     0.00     7.59     1.57
3:     0.00     0.31     0.15     0.00     0.35     0.04     0.00     0.73     0.13
4:     5.20     0.13     5.07     3.73     1.29     2.47     9.82     1.57     8.32
5:     0.22     0.07     0.00     0.17     0.02     0.00     0.22     0.09     0.00
6:     0.00     7.99     3.06     0.00     5.16     0.24     0.00     2.90     0.71
7:     0.00     0.34     0.07     0.00     0.51     0.02     0.00     0.92     0.09
8:    11.24     3.06     8.33     5.76     0.24     5.68     4.40     0.71     3.82
9:     0.24     0.02     0.00     0.06     0.08     0.00     0.07     0.13     0.00
10:     0.00     4.54     1.77     0.00     3.40     3.24     0.00     7.43     0.42
11:     0.00     0.38     0.02     0.00     0.85     0.08     0.00     0.06     0.13
12:     6.48     1.77     4.92     7.50     3.24     4.25     7.98     0.42     7.49
13:     0.03     0.12     0.00     0.05     0.03     0.00     0.13     0.12     0.00
14:     0.00     4.61     2.49     0.00     1.40     3.20     0.00     7.20     1.06
15:     0.00     0.59     0.12     0.00     0.79     0.03     0.00     0.53     0.12
16:     7.78     2.49     5.19     5.37     3.20     2.19     8.78     1.06     7.73

M lhs matrix:
1:     1.46     0.89     0.30     1.30     0.77     0.03     0.21
2:     0.89    83.17     0.00   168.15     1.26     0.06     0.32
3:     0.30     0.00     1.81     3.71     1.07     0.05     0.30
4:     1.30   168.15     3.71   413.88     0.68     0.03     0.13
5:     0.77     1.26     1.07     0.68     0.90     0.00     0.00
6:     0.03     0.06     0.05     0.03     0.00     0.90     0.00
7:     0.21     0.32     0.30     0.13     0.00     0.00     0.90
M rhs matrix:
1:     0.26     0.27     0.26     0.13
2:     5.14     8.85     5.77     6.57
3:     0.48     0.44     0.41     0.66
4:    13.71    21.53    15.63    17.63
5:     0.91     0.31     0.97     0.56
6:     0.17     0.38     0.15     0.57
7:     0.57     0.14     0.60     0.80
P lhs matrix:
1:   204.92   144.94   177.23     0.98     1.14     0.74     0.80
2:   144.94   161.29   117.57     0.67     0.76     0.52     0.58
3:   177.23   117.57   184.56     0.78     0.99     0.68     0.66
4:     0.98     0.67     0.78     0.90     0.00     0.00     0.00
5:     1.14     0.76     0.99     0.00     0.90     0.00     0.00
6:     0.74     0.52     0.68     0.00     0.00     0.90     0.00
7:     0.80     0.58     0.66     0.00     0.00     0.00     0.90
P rhs matrix:
1:    16.50    11.64    19.23
2:    15.04    11.92    11.75
3:    16.01     8.79    18.08
4:     0.91     0.17     0.57
5:     0.31     0.38     0.14
6:     0.97     0.15     0.60
7:     0.56     0.57     0.80

*/
int main ( void){
    const uint32_t ncluster = 5;
    const uint32_t ncycle = 3;
    MAT matWe = coerce_MAT_from_array(ncluster,1,we);
    fputs("Weight matrix:\n",stdout);
    show_MAT(xstdout,matWe,0,0);
    MAT matLambda = coerce_MAT_from_array(ncluster,1,lambda);
    fputs("Lambda matrix:\n",stdout);
    show_MAT(xstdout,matLambda,0,0);
    MAT matInts = coerce_MAT_from_array(NBASE*ncycle,ncluster,ints_t);
    fputs("Intensity matrix:\n",stdout);
    show_MAT(xstdout,matInts,0,0);
    
    MAT J = calculateJ(matLambda,matWe,bases,ncycle,NULL);
    MAT Jt = transpose(J);
    fputs("J matrix:\n",stdout);
    show_MAT(xstdout,J,0,0);
    MAT K = calculateK(matLambda,matWe,bases,matInts,ncycle,NULL);
    MAT Kt = transpose(K);
    fputs("K matrix:\n",stdout);
    show_MAT(xstdout,K,0,0);
    
    real_t * tmp = calloc(NBASE*NBASE+ncycle*ncycle,sizeof(real_t));
    MAT matSbarT = coerce_MAT_from_array(ncycle,NBASE,SbarT);
    MAT matSbar = transpose(matSbarT);
    fputs("SbarT matrix:\n",stdout);
    show_MAT(xstdout,matSbarT,0,0);
    MAT matM = coerce_MAT_from_array(NBASE,NBASE,M);
    MAT matMt = transpose(matM);
    fputs("M matrix:\n",stdout);
    show_MAT(xstdout,matM,0,0);
    MAT matP = coerce_MAT_from_array(ncycle,ncycle,P);
    fputs("P matrix:\n",stdout);
    show_MAT(xstdout,matP,0,0);

    MAT matIbar = coerce_MAT_from_array(NBASE,ncycle,Ibar);
    MAT matIbarT = transpose(matIbar);
    fputs("Ibar matrix:\n",stdout);
    show_MAT(xstdout,matIbar,0,0);

    fputs("M lhs matrix:\n",stdout);
    MAT lhs = calculateMlhs(0.9,matSbarT,matP,Jt,tmp,NULL);
    show_MAT(xstdout,lhs,0,0);
    MAT rhs = calculateMrhs(matIbarT,matP,Kt,tmp,NULL);
    fputs("M rhs matrix:\n",stdout);
    show_MAT(xstdout,rhs,0,0);
    
    int retM = solver(lhs,rhs);
    xfprintf(xstdout,"M solution matrix (info=%d):\n",retM);
    show_MAT(xstdout,rhs,0,0);

    MAT Plhs = calculatePlhs(0.9,matSbar,matMt,J,tmp,NULL);
    fputs("P lhs matrix:\n",stdout);
    show_MAT(xstdout,Plhs,0,0);
    MAT Prhs = calculatePrhs(matIbar,matMt,K,tmp,NULL);
    fputs("P rhs matrix:\n",stdout);
    show_MAT(xstdout,Prhs,0,0);

    int retP = solver(lhs,rhs);
    xfprintf(xstdout,"P solution matrix (info=%d):\n",retM);
    show_MAT(xstdout,rhs,0,0);
    
    return EXIT_SUCCESS;
}
#endif

