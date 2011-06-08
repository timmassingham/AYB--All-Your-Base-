/*
 *  Copyright (C) 2010 by Tim Massingham, European Bioinformatics Institute
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

/* Standard copyright header */
/* Description
 * Include statement about row/col major format
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <err.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "matrix.h"
#include "lapack.h"

#define WARN_MEM(A) warn("Failed to allocation memory for %s at %s:%d.\n",(A),__FILE__,__LINE__)

/*  First deal with allocation and deallocation of matrices  */
/* Allocate memory for matrix of a specified size */
MAT new_MAT( const int nrow, const int ncol ){
    MAT mat = malloc(sizeof(*mat));
    if ( NULL==mat){
        WARN_MEM("matrix");
        return NULL;
    }
    mat->ncol=ncol;
    mat->nrow=nrow;
    /* Number of rows or columns might be zero but probably an error, so warn
     * as such. Want to avoid malloc(0) since this is "implementation defined"
     * in the C standard, may be a real address that should not be used or NULL
     */
     if ( 0==ncol || 0==nrow ){
         warn("One of dimensions of matrix equal to zero. Setting memory to NULL.\n");
         mat->x = NULL;
     } else {
         /* Usual case. Use calloc rather than malloc so elements are
          * initialised
          */
         mat->x = calloc(nrow*ncol,sizeof(real_t));
         if ( NULL==mat->x){
             WARN_MEM("matrix elements");
             free(mat);
             mat = NULL;
         }
     }
     
     return mat;
}

MAT new_MAT_from_array( const uint32_t nrow, const uint32_t ncol, const real_t * x){
    if(NULL==x){ return NULL;}
    MAT mat = new_MAT(nrow,ncol);
    if(NULL==mat){return NULL;}
    memcpy(mat->x,x,nrow*ncol*sizeof(real_t));
    return mat;
}

/* Free memory allocated for matrix */
void free_MAT ( MAT mat ){
    if(NULL==mat){ return; }
    /* Memory for elements may be NULL if nrow or ncol equals zero */
    if ( NULL!=mat->x){
        free(mat->x);
    }
    free(mat);
}

// Frees matrix but not elements
void shallow_free_MAT( MAT mat){
	if(NULL!=mat){ free(mat);}
}

MAT coerce_MAT_from_array(const uint32_t nrow, const uint32_t ncol, real_t * x){
	assert(NULL!=x);
	MAT mat = malloc(sizeof(*mat));
	if(NULL==mat){ return NULL; }
	mat->nrow = nrow;
	mat->ncol = ncol;
	mat->x = x;
	return mat;
}


void show_MAT ( XFILE * fp, const MAT mat, const uint32_t mrow, const uint32_t mcol){
    if(NULL==fp){ return;}
    if(NULL==mat){ return;}
    
    const uint32_t nrow = mat->nrow;
    const uint32_t ncol = mat->ncol;
    const uint32_t maxrow = (mrow!=0 && mrow<nrow)?mrow:nrow;
    const uint32_t maxcol = (mcol!=0 && mcol<ncol)?mcol:ncol;
    xfprintf(fp,"%d\t%d\n",maxrow,maxcol);
    for( int row=0 ; row<maxrow ; row++){
        for ( int col=0 ; col<maxcol ; col++){
            xfprintf(fp," %#8.6e",mat->x[col*nrow+row]);
        }
        if(maxcol<ncol){ xfprintf(fp,"\t... (%u others)",ncol-maxcol); }
        xfputc('\n',fp);
    }
    if( maxrow<nrow){ xfprintf(fp,"... (%u others)\n",nrow-maxrow); }
}

MAT copy_MAT( const MAT mat){
    if(NULL==mat){ return NULL;}

    MAT newmat = new_MAT(mat->nrow,mat->ncol);
    if(NULL==newmat){ return NULL;}
    
    memcpy(newmat->x,mat->x,mat->nrow*mat->ncol*sizeof(real_t));
    return newmat;
}

MAT copyinto_MAT( MAT newmat, const MAT mat){
    if(newmat->nrow!=mat->nrow){ return NULL;}
    if(newmat->ncol!=mat->ncol){ return NULL;}
    memcpy(newmat->x,mat->x,mat->nrow*mat->ncol*sizeof(real_t));
    return newmat;
}

MAT set_MAT( MAT mat, const real_t x){
	if(NULL==mat){ return NULL;}
	const uint32_t nelt = mat->nrow * mat->ncol;
	for ( uint32_t i=0 ; i<nelt ; i++){
		mat->x[i] = x;
	}
	return mat;
}


bool is_square(const MAT mat){
    validate(NULL!=mat,false);
    if(mat->nrow!=mat->ncol){ return false;}
    return true;
}


// Copy lower diagonal of matrix to upper diagonal
MAT symmeteriseL2U( MAT mat){
    validate(NULL!=mat,NULL);
    validate(mat->nrow==mat->ncol,NULL);
    const uint32_t n = mat->ncol;
    for ( uint32_t col=0 ; col<n ; col++){
        for ( uint32_t row=col ; row<n ; row++){
            mat->x[row*n+col] = mat->x[col*n+row];
        }
    }
    return mat;
}


MAT identity_MAT( const int nrow){
    MAT mat = new_MAT(nrow,nrow);
    validate(NULL!=mat,NULL);
    for ( int i=0 ; i<nrow ; i++){
        mat->x[i*nrow+i] = 1.0;
    }
    return mat;
}

MAT reshape_MAT( MAT mat, const int nrow){
    validate(NULL!=mat,NULL);
    validate(((mat->nrow*mat->ncol)%nrow)==0,NULL);
    mat->ncol = (mat->nrow*mat->ncol)/nrow;
    mat->nrow = nrow;
    return mat;
}

MAT trim_MAT( MAT mat, const int mrow, const int mcol, const bool forwards){
    validate(NULL!=mat,NULL);
    validate(mrow>=0,NULL);
    validate(mcol>=0,NULL);
    validate(mrow<=mat->nrow,NULL);
    validate(mcol<=mat->ncol,NULL);
    if(forwards==false){ errx(EXIT_FAILURE,"Forwards==false not implemented in %s (%s:%d)\n",__func__,__FILE__,__LINE__);}
    for ( uint32_t col=0 ; col<mcol ; col++){
        uint32_t midx = col*mrow;
        uint32_t nidx = col*mat->nrow;
        memmove(mat->x+midx,mat->x+nidx,mrow*sizeof(real_t));
    }
    mat->nrow = mrow;
    mat->ncol = mcol;
    return mat;
}

MAT * block_diagonal_MAT( const MAT mat, const int n){
    validate(NULL!=mat,NULL);
    validate(mat->ncol==mat->nrow,NULL);    // Ensure symmetry
    const int nelts = mat->ncol / n;        // Number of blocks on diagonal
    validate((mat->ncol % n)==0,NULL);      // Is parameter valid?
    
    // Create memory
    MAT * mats = calloc(nelts,sizeof(*mats));
    if(NULL==mats){ goto cleanup; }
    for ( uint32_t i=0 ; i<nelts ; i++){
        mats[i] = new_MAT(n,n);
        if(NULL==mats[i]){ goto cleanup;}
    }
    // Copy into diagonals
    for ( uint32_t i=0 ; i<nelts ; i++){
        for ( uint32_t col=0 ; col<n ; col++){
            const uint32_t oldcol = i*n+col;
            for ( uint32_t row=0 ; row<n ; row++){
                const uint32_t oldrow = i*n+row;
                mats[i]->x[col*n+row] = mat->x[oldcol*mat->nrow+oldrow];
            }
        }
    }
    return mats;
    
cleanup:
    if(NULL!=mats){
        for ( uint32_t i=0 ; i<nelts ; i++){
            free_MAT(mats[i]);
        }
    }
    xfree(mats);
    return NULL;
}

MAT scale_MAT(MAT mat, const real_t f){
    validate(NULL!=mat,NULL);
    const uint32_t nelt = mat->ncol * mat->nrow;
    for ( uint32_t elt=0 ; elt<nelt ; elt++){
            mat->x[elt] *= f;
    }
    return mat;
}

// Only works for square matrices
MAT transpose_inplace( MAT mat){
    validate(NULL!=mat,NULL);
    const uint32_t ncol = mat->ncol;
    const uint32_t nrow = mat->nrow;
    validate(ncol==nrow,NULL);
    
    for ( uint32_t col=0 ; col<ncol ; col++){
        for ( uint32_t row=0 ; row<col ; row++){
            real_t x = mat->x[col*nrow + row];
            mat->x[col*nrow + row] = mat->x[row*ncol + col];
            mat->x[row*ncol + col] = x;
        }
    }
    mat->nrow = ncol;
    mat->ncol = nrow;
    return mat;
}

MAT transpose( const MAT mat){
    validate(NULL!=mat,NULL);
    MAT tmat = new_MAT(mat->ncol,mat->nrow);
    validate(NULL!=tmat,NULL);
    for ( uint32_t col=0 ; col<mat->ncol ; col++){
        for ( uint32_t row=0 ; row<mat->nrow ; row++){
            tmat->x[row*mat->ncol+col] = mat->x[col*mat->nrow+row];
        }
    }
    return tmat;
}

/* Invert a symmetric matrix
 */
MAT invertSym(const MAT mat){
   validate(NULL!=mat,NULL);
   validate(mat->nrow==mat->ncol,NULL);

   MAT A = copy_MAT(mat);
   if(NULL==A){ return NULL;}
   int N = mat->nrow;
   int INFO;
   // Cholesky
   potrf(LAPACK_UPPER,&N,A->x,&N,&INFO);
   instrument(fprintf(stderr,"potrf in %s returned %d\n",__func__,INFO));
   // Invert
   potri(LAPACK_UPPER,&N,A->x,&N,&INFO);
   instrument(fprintf(stderr,"potri in %s returned %d\n",__func__,INFO));
   // Make Symmetric
   for ( int i=0 ; i<N ; i++){
      for (int j=0 ; j<i ; j++){
	 A->x[j*N+i] = A->x[i*N+j];
      }
   }
   return A;
}


MAT invert(const MAT mat){
    validate(NULL!=mat,NULL);
    validate(mat->nrow==mat->ncol,NULL);
    
    int INFO;
    int N;
    int * IPIV=NULL;
    int LWORK=0;
    double * WORK=NULL,WORKSIZE=0;
    
    N = mat->nrow;
    // Get temporary memory for inversion
    LWORK = -1;
    getri(&N,mat->x,&N,IPIV,&WORKSIZE,&LWORK,&INFO);
    LWORK=(int)WORKSIZE;
    WORK = malloc(LWORK*sizeof(double));
    IPIV = malloc(N*sizeof(int));
    
    MAT matinv = copy_MAT(mat);
    // LU decomposition required for inversion
    getrf(&N,&N,matinv->x,&N,IPIV,&INFO);
    instrument(fprintf(stderr,"getrf in %s returned %d\n",__func__,INFO));
    // Invert
    getri(&N,matinv->x,&N,IPIV,WORK,&LWORK,&INFO);
    instrument(fprintf(stderr,"getri in %s returned %d\n",__func__,INFO));
    
    free(IPIV);
    free(WORK);
    
    return matinv;
}

real_t xMy( const real_t * x, const MAT M, const real_t * y){
    validate(NULL!=x,NAN);
    validate(NULL!=M,NAN);
    validate(NULL!=y,NAN);
    const uint32_t ncol = M->ncol;
    const uint32_t nrow = M->nrow;
    
    real_t res = 0.;
    for ( uint32_t col=0 ; col<ncol ; col++){
        real_t rowtot = 0.;
        for ( uint32_t row=0 ; row<nrow ; row++){
            rowtot += x[row] * M->x[col*nrow+row];
        }
        res += rowtot * y[col];
    }
    return res;
}


real_t normalise_MAT(MAT mat, const real_t delta_diag){
    validate(NULL!=mat,NAN);
    validate(mat->nrow==mat->ncol,NAN);
    const int n = mat->nrow;
    int piv[n];
    
    if(0.0!=delta_diag){
        for ( uint32_t i=0 ; i<n ; i++){
            mat->x[i*n+i] += delta_diag;
        }
    }
    
    MAT mcopy = copy_MAT(mat);
    int info = 0;
    getrf(&n,&n,mcopy->x,&n,piv,&info);
    instrument(fprintf(stderr,"getrf in %s returned %d\n",__func__,info));
    
    real_t logdet = 0.;
    for ( uint32_t i=0 ; i<n ; i++){
        logdet += log(fabs(mcopy->x[i*n+i]));
    }
    free_MAT(mcopy);
    
    real_t f = 1e-5 + exp(logdet/n);
    scale_MAT(mat,1./f);
    return f;
}

MAT cholesky( MAT mat){
    validate(NULL!=mat,NULL);
    validate(mat->nrow==mat->ncol,NULL);
    int info=0;
    int n = mat->nrow;
    potrf(LAPACK_UPPER,&mat->nrow,mat->x,&mat->nrow,&info);
    instrument(fprintf(stderr,"potrf in %s returned %d\n",__func__,info));
    for( int i=0 ; i<n ; i++){
        for ( int j=i+1 ;j<n ; j++){
           mat->x[i*n+j] = 0;
        }
    } 
    return mat;
}

struct structLU LUdecomposition( const MAT mat){
   struct structLU structLUnill = {NULL,NULL};
   validate(NULL!=mat,structLUnill);
   validate(mat->nrow==mat->ncol,structLUnill);
   const int n = mat->nrow;
   int * piv = calloc(n,sizeof(int));;

   MAT mcopy = copy_MAT(mat);
   int info = 0;
   getrf(&n,&n,mcopy->x,&n,piv,&info);

   return (struct structLU) { mcopy, piv};
}

