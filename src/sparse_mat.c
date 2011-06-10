#include <math.h>
#include <strings.h>
#include "sparse_mat.h"

SparseMAT make_SparseMAT(const MAT m, const real_t thresh){
	if(NULL==m){ return NULL; }
	if(isnan(thresh) || thresh<0){ return NULL;}
	
	// Count elements
	int nelt = 0;
	int celt[m->ncol];
	bzero(celt,m->ncol*sizeof(int));
	for ( int col=0,idx=0 ; col<m->ncol ; col++){
		for ( int row=0 ; row<m->nrow ; row++,idx++){
			if(fabs(m->x[idx])>thresh){
				nelt++;
				celt[col]++;
			}
		}
	}

	// Create structure
	SparseMAT spm = calloc(1,sizeof(*spm));
	if(NULL==spm){ return NULL;}
	spm->rowidx = calloc(nelt+m->ncol,sizeof(int));
	if(NULL==spm->rowidx){ goto cleanup; }
	spm->val = calloc(nelt,sizeof(real_t));
	if(NULL==spm->val){ goto cleanup; }
	spm->nrow = m->nrow;
	spm->ncol = m->ncol;
	spm->nelt = nelt;

	// Fill structure
	for ( int col=0,idx=0,i=0 ; col<m->ncol ; col++){
		spm->rowidx[col+idx] = celt[col];
		for ( int row=0 ; row<m->nrow ; row++,i++){
			if(fabs(m->x[i])>thresh){
				spm->val[idx] = m->x[i];
				spm->rowidx[col+idx+1] = row;
				idx++;
			}
		}
	}
	fprintf(stderr,"Number of sparse elements = %d / %d = %f%%\n",spm->nelt,m->nrow*m->ncol,(100.0*spm->nelt)/(m->nrow*m->ncol));

	return spm;

cleanup:
	free_SparseMAT(spm);
	return NULL;
}

void free_SparseMAT( SparseMAT spm){
	if(NULL==spm){ return; }
	free(spm->rowidx);
	free(spm->val);
	free(spm);
}

real_t * sparseMv( const SparseMAT spm, const real_t * v, real_t * res){
	if(NULL==spm){ return NULL; }
	if(NULL==v){ return NULL; }
	if(NULL==res){
		res = calloc(spm->nrow,sizeof(real_t));
		if(NULL==res){ return NULL; }
	}
	bzero(res,spm->nrow*sizeof(real_t));

	for ( int col=0,idx=0 ; col<spm->ncol ; col++){
		int ncelt = spm->rowidx[col+idx];
		for ( int celt=0 ; celt<ncelt ; celt++,idx++){
			int row = spm->rowidx[col+idx+1];
			res[row] += spm->val[idx] * v[col];
		}
	}

	return res;
}
