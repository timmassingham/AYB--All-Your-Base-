#ifndef SPARSE_MAT_H
#define SPARSE_MAT_H

#include "matrix.h"
#include "utility.h"

typedef struct {
                int nrow,ncol,nelt;
                int * rowidx;
                real_t * val;
} * SparseMAT;

SparseMAT make_SparseMAT(const MAT m, const real_t thresh);
void free_SparseMAT( SparseMAT spm);
real_t * sparseMv( const SparseMAT spm, const real_t * v, real_t * res);


#endif /* SPARSE_MAT_H */

