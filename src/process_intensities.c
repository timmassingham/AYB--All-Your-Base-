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

#include "matrix.h"
#include "sparse_mat.h"
#include "nuc.h"
#include <strings.h>
#include <math.h>
#include "lapack.h"


real_t calculateLSSi_sub( const real_t lambda, const NUC * bases, const MAT M, const MAT P, const MAT N, const real_t * I, MAT e);
MAT expected_intensities( const real_t lambda, const NUC * bases, const MAT M, const MAT P, const MAT N, MAT e);

/*
 * Process intensities, p = Minv %*% (Intensities-N) %*% Pinv
 * Uses identity: Vec(p) = ( Pinv^t kronecker Minv) Vec(Intensities-N)
 * Storing Intensities-N as an intermediate saved < 3%
 * Calculating p^t rather than p (pcol loop is over minor index) made no difference
 * Using Pinv rather than Pinv^t makes little appreciable difference.
 */
MAT process_intensities(const int16_t * intensities, const MAT Minv_t, const MAT Pinv_t, const MAT N, MAT p){
    validate(NULL!=intensities,NULL);
    validate(NULL!=Minv_t,NULL);
    validate(NULL!=Pinv_t,NULL);
    validate(NULL!=N,NULL);
    
    const uint_fast32_t ncycle = Pinv_t->nrow;
    if(NULL==p){
        p = new_MAT(NBASE,ncycle);
        validate(NULL!=p,NULL);
    }
    bzero(p->x,p->nrow * p->ncol * sizeof(real_t));
    
    for( uint_fast32_t icol=0 ; icol<ncycle ; icol++){    // Columns of Intensity
        real_t dp[NBASE] = {0,0,0,0};
        for( uint_fast32_t base=0 ; base<NBASE ; base++){ // Bases (rows of Minv, cols of Minv_t)
            for ( uint_fast32_t chan=0 ; chan<NBASE ; chan++){  // Channels
                dp[base] += Minv_t->x[base*NBASE+chan] * (intensities[icol*NBASE+chan] - N->x[icol*NBASE+chan]);
            }
        }
        for ( uint_fast32_t pcol=0 ; pcol<ncycle ; pcol++){ // Columns of p
            const real_t tmp = Pinv_t->x[icol*ncycle+pcol];
            for( uint_fast32_t base=0 ; base<NBASE ; base++){
                p->x[pcol*NBASE+base] += tmp * dp[base];
            }
        }
    }

    return p;
}

MAT processNew( const struct structLU AtLU, const MAT N, const int16_t * intensities, MAT p){
	if(NULL==AtLU.mat || NULL==N || NULL==intensities){ return NULL;}
	const int ncycle = N->ncol;
	if(NULL==p){
		p = new_MAT(4,ncycle);
		if(NULL==p){ return NULL;}
	}
	const int nelt = 4*ncycle;

	//real_t intminusN[nelt];
	for ( int i=0 ; i<nelt ; i++){
	   p->x[i] = intensities[i]-N->x[i];
	}

	const int inc = 1;
	int info = 0;
	getrs(LAPACK_TRANS,&nelt,&inc,AtLU.mat->x,&nelt,AtLU.piv,p->x,&nelt,&info);
	//instrument(fprintf(stderr,"getrs returned %d in %s\n",info,__func__));
	return p;
}

MAT processSparseNew( const SparseMAT spAinv, const MAT N, const int16_t * intensities, MAT p){
	if(NULL==spAinv || NULL==N || NULL==intensities){ return NULL;}
	const int ncycle = N->ncol;
	if(NULL==p){
		p = new_MAT(4,ncycle);
		if(NULL==p){ return NULL;}
	}

	const int nelt = 4*ncycle;
	real_t intminusN[nelt];
	for ( int i=0 ; i<nelt ; i++){
		intminusN[i] = intensities[i]-N->x[i];
	}

	sparseMv(spAinv,intminusN,p->x);
	return p;
}


MAT expectedNew(const MAT A, const MAT N, const NUC * bases, MAT e){
	if( NULL==A || NULL==N || NULL==bases){ return NULL; }
	const int ncycle = N->ncol;
	const int lda = 4*ncycle;
	if(NULL==e){
		e = new_MAT(4,ncycle);
		if(NULL==e){ return NULL; }
	}
	bzero(e->x,lda*sizeof(real_t));

	// A vec(S)
	for ( int cy=0 ; cy<ncycle ; cy++){
		const int idx = cy*4 + bases[cy];
		for ( int i=0 ; i<lda ; i++){
			e->x[i] += A->x[idx*lda+i];
		}
	}

	// + N
	for ( int i=0 ; i<lda ; i++){
		e->x[i] += N->x[i];
	}

	return e;
}

real_t calculateLSS( const MAT lambda, const NUC * bases, const MAT M, const MAT P, const MAT N, const MAT I){
    validate(NULL!=lambda,NAN);
    validate(NULL!=bases,NAN);
    validate(NULL!=M,NAN);
    validate(NULL!=P,NAN);
    validate(NULL!=N,NAN);
    validate(NULL!=I,NAN);
    const uint32_t ncycle = P->nrow;
    const uint32_t ncluster = I->ncol;
    
    MAT e = new_MAT(NBASE,ncycle);
    validate(NULL!=e,NAN);
    real_t lss = 0;
    for ( uint32_t i=0 ; i<ncluster ; i++){
        lss += calculateLSSi_sub(lambda->x[i],bases+i*ncycle,M,P,N,I->x+i*ncycle*NBASE,e);
    }
    free_MAT(e);
    return lss;
}

real_t calculateWeLSS( const MAT we, const MAT lambda, const NUC * bases, const MAT M, const MAT P, const MAT N, const MAT I){
    validate(NULL!=we,NAN);
    validate(NULL!=lambda,NAN);
    validate(NULL!=bases,NAN);
    validate(NULL!=M,NAN);
    validate(NULL!=P,NAN);
    validate(NULL!=N,NAN);
    validate(NULL!=I,NAN);
    const uint32_t ncycle = P->nrow;
    const uint32_t ncluster = I->ncol;
    
    MAT e = new_MAT(NBASE,ncycle);
    validate(NULL!=e,NAN);
    real_t lss = 0;
    for ( uint32_t i=0 ; i<ncluster ; i++){
        lss += we->x[i] * calculateLSSi_sub(lambda->x[i],bases+i*ncycle,M,P,N,I->x+i*ncycle*NBASE,e);
    }
    free_MAT(e);
    return lss;
}

real_t calculateLSSi_sub( const real_t lambda, const NUC * bases, const MAT M, const MAT P, const MAT N, const real_t * I, MAT e){
    validate(lambda>=0,NAN);
    validate(NULL!=bases,NAN);
    validate(NULL!=M,NAN);
    validate(NULL!=P,NAN);
    validate(NULL!=N,NAN);
    validate(NULL!=I,NAN);
    validate(NULL!=e,NAN);
    const uint32_t ncycle = P->nrow;
    
    expected_intensities(lambda,bases,M,P,N,e);
    real_t lss = 0;
    for ( uint32_t i=0 ; i<(ncycle*NBASE) ; i++){
        real_t tmp = (I[i] - e->x[i]);
        lss += tmp*tmp;
    }
    return lss;
}

MAT expected_intensities( const real_t lambda, const NUC * bases, const MAT M, const MAT P, const MAT N, MAT e){
    validate(lambda>=0,NULL);
    validate(NULL!=bases,NULL);
    validate(NULL!=M,NULL);
    validate(NULL!=P,NULL);
    validate(NULL!=N,NULL);
    const uint_fast32_t ncycle = P->nrow;
    if(NULL==e){
        e = new_MAT(NBASE,ncycle);
        validate(NULL!=e,NULL);
    }
    bzero(e->x,NBASE*ncycle*sizeof(real_t));
    for(uint_fast32_t cy2=0 ; cy2<ncycle ; cy2++){
        for(uint_fast32_t cy=0 ; cy<ncycle ; cy++){
            const int base = bases[cy];
            for ( uint_fast32_t ch=0 ; ch<NBASE ; ch++){
                e->x[cy2*NBASE+ch] += M->x[base*NBASE+ch] * P->x[cy2*ncycle+cy];
            }
        }
    }
        
    // Multiply by brightness;
    scale_MAT(e,lambda);
    // Add noise
    for ( uint_fast32_t i=0 ; i<(NBASE*ncycle) ; i++){
        e->x[i] += N->x[i];
    }
    return e;
}
        

#ifdef TEST
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "xio.h"

const uint32_t ncycle = 5;
const real_t sMinv_t[] = {
    0.8901800, 0.2411217, 0.6736770, 0.6388693,
    0.8930780, 0.6432912, 0.1419963, 0.6806006,
    0.3074286, 0.6463750, 0.5734479, 0.1899474,
    0.1081050, 0.2802771, 0.1109772, 0.1573850
};

const real_t sI[] = {
    0.7800604, 0.0900520, 0.6934569, 0.9832249,
    0.6007023, 0.7111706, 0.3739851, 0.7511164,
    0.1806179, 0.3007264, 0.5049830, 0.5060436,
    0.5215840, 0.4716412, 0.3251011, 0.8108462,
    0.1337858, 0.5297113, 0.6538163, 0.7541822
};

const real_t sN[] = {
    0.9297764, 0.07949248, 0.09707808, 0.67001314,
    0.6692831, 0.26673910, 0.18551665, 0.47152955,
    0.7440895, 0.75429913, 0.39306556, 0.54502925,
    0.9999948, 0.52553739, 0.44663278, 0.62640211,
    0.5235454, 0.30976778, 0.31814943, 0.08174834
};

const real_t sPinv_t[] = {
    0.587633975, 0.2880588, 0.9676281, 0.1641735, 0.5512092,
    0.544311555, 0.6446207, 0.7931918, 0.6183713, 0.9766380,
    0.002392489, 0.5103488, 0.4066239, 0.8455400, 0.3015657,
    0.218021307, 0.8041612, 0.7649060, 0.5440081, 0.7128468,
    0.219031750, 0.1337430, 0.8726016, 0.9766849, 0.9514125
};

/*  Result should be:
 * 0.4583560 -0.19921692 0.5144813 -0.04488226 0.4911784
 * 0.3272714 -0.32163246 0.1783205 -0.28063772 0.3147591
 * 0.4723463  0.04245129 0.6562353  0.19414411 0.6650795
 * 0.1831871  0.03456732 0.2721046  0.10659499 0.2975027
 */
    
int main (int argc, char * argv[]){
    MAT Minv_t, Pinv_t,I,N;

    if(2!=argc){ fputs("Usage: test ntimes\n",stdout); return EXIT_FAILURE; }
    unsigned int ntimes;
    sscanf(argv[1],"%u",&ntimes);

    Minv_t = new_MAT(NBASE,NBASE);
    memcpy(Minv_t->x,sMinv_t,NBASE*NBASE*sizeof(real_t));
    Pinv_t = new_MAT(ncycle,ncycle);
    memcpy(Pinv_t->x,sPinv_t,ncycle*ncycle*sizeof(real_t));
    I = new_MAT(NBASE,ncycle);
    memcpy(I->x,sI,NBASE*ncycle*sizeof(real_t));
    N = new_MAT(NBASE,ncycle);
    memcpy(N->x,sN,NBASE*ncycle*sizeof(real_t));
    
    
    fputs("Minv_t:\n",stdout);
    show_MAT(xstdout,Minv_t,NBASE,NBASE);
    fputs("Pinv:\n",stdout);
    show_MAT(xstdout,Pinv_t,ncycle,ncycle);
    fputs("N:\n",stdout);
    show_MAT(xstdout,N,NBASE,ncycle);
    fputs("I:\n",stdout);
    show_MAT(xstdout,I,NBASE,ncycle);
    
    
    fputs("\nResult:\n",stdout);
    MAT proc = NULL;
    for ( uint32_t i=0 ; i<ntimes ; i++){
        proc = process_intensities(I->x,Minv_t,Pinv_t,N,proc);
    }
    show_MAT(xstdout,proc,NBASE,ncycle);
    
    return EXIT_SUCCESS;
};

#endif
