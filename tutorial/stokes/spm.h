#ifndef _SPM_H_
#define _SPM_H_

typedef void* spm;

extern double* spm_double(spm *A, long i, long j);
extern long*   nonzero_col(spm *p);

static long* nonzero_n;
static long  nonzero_k;
#define fornonzero(Ai,j) \
for(nonzero_n=nonzero_col(Ai),nonzero_k=dim1(nonzero_n), \
j=nonzero_n[nonzero_k];nonzero_k;j=nonzero_n[--nonzero_k])

#endif
