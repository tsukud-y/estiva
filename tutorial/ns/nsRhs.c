#include <stdio.h>

#include "estiva/ary.h"

#include "ns.h"
#include "fem.h"

void foo_mulmx(double **tp, MX *A, double *x){
  long i, j, J, m = A->m, n = A->n;
  double *t;
  
  ary1(*tp, m+1);
  
  t  = *tp;
  
  for(i=0; i<=m; i++) t[i] = 0.0;


  
  /*
  for( i = 1; i <= m; i++ ) for ( j = 1; j <= m; j++ ) {
      t[  i] += mx(A,i,j)*x[  j];
    }
  */

  mx(A,1,1) = mx(A,1,1);

  for(i=0; i< m; i++) for(j=0; j< n; j++) {
      J = A->IA[i][j];
      if (J) t[i+1] += A->A[i][j]*x[J];
    }
}



void estiva_nsRhs(double *b, xyc *Z, nde *N, MX *M, double *x, double t)
{
  static count = 1;
  long   i, j, NUM, m, n;


  m = dimp2(N);
  n = dim1(Z);

  NUM = m*2+n;

  for( i = 1; i <= NUM; i++ ) b[i] = 0.0;



  for( i = 1; i <= m; i++ ) for ( j = 1; j <= m; j++ ) {
      b[  i] += mx(M,i,j)*x[  j];
      b[m+i] += mx(M,i,j)*x[m+j];
    }
#if 0
  {
    static int init = 0;
    static MX *M2;
    if ( !init ) {
      long i, j, m, n;
      init = 1;
      m = dimp2(N);
      n = dim1(Z);
      initmx(M2,NUM,50);
      for ( i = 1; i <= m; i++ ) for ( j = 1; j <=m; j++ ) {
	  mx(M2,  i,  j) = mx(M,i,j);
	  mx(M2,m+i,m+j) = mx(M,i,j);
	}
      pltmx(M2);pltmx(M);
    }
    mulmx(b,M2,x);
  }
#endif
  {
    static double *bx, *by, *xx, *xy;
    ary1(xx,m+1); ary1(xy,m+1);
    for ( i = 1; i<= m; i++ ) {
      xx[i] = x[i];
      xy[i] = x[i+m];
    }
    foo_mulmx(&bx,M,xx);
    foo_mulmx(&by,M,xy);
    //mulmx(bx,M,xx);
    //mulmx(by,M,xy);
    for(i=1;i<=m;i++){
      b[i] = bx[i];
      b[i+m] = by[i];
    }
  }
}
