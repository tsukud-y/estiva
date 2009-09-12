#include <stdio.h>
#include "ary.h"
#include "spm.h"

#define max(x,y) (x>y?x:y)
#define A(i,j) mx(A,i,j)


static long count_m(void *A)
{
  long i, j, n, m;

  n = dim1(A);
  m = 0;

  for(i=1;i<=n;i++) for(j=1;j<=n;j++) if(A(i,j) != 0.0) m = max(m,j);

  return m;
}

void cplot(void* A)
{
  long i,j,n,m;
  n = dim1(A);
  m = count_m(A);

  for(i=1;i<=n;i++){
    for(j=1;j<=m;j++)
      if(A(i,j) == 1.0) printf("£±");
      else if(A(i,j) != 0.0) printf("¢£");
      else  printf("¢¢");
    printf("\n");
  }
}
