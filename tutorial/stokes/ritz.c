#include <stdio.h>
#include "ary.h"
#include "spm.h"


#define A(i,j) mx(A,i,j)


static int this_is_boundary(spm *A, long j)
{
  int  i;
  long n;

  n = dim1(A);

  if(A(j,j) != 1.0) return 0;

  for(i=1;i<=n;i++) if(A(j,i) != 0.0 && j != i) return 0;
  return 1;
}


void ritz(spm *A, double *b)
{
  static long *P;
  int i, j, n, k;
  n = dim1(A);

  ary1(P,n+1);
  k=0;
  for(i=1;i<=n;i++){
    if(this_is_boundary(A,i)) P[k++] = i; 
  }
  P[k] = 0;

  for(k=0;P[k];k++){
    j = P[k];
    for(i=1;i<=n;i++) A(i,j) *= b[j];
    A(j,j) = 0.0;
  }
  for(k=0;P[k];k++){
    j = P[k];
    for(i=1;i<=n;i++) b[i] -= A(i,j);
  }
  for(k=0;P[k];k++){
    j = P[k];
    for(i=1;i<=n;i++) A(i,j) = 0.0;
    A(j,j) = 1.0;
  }
}
