#include "ns.h"
#include "fem.h"

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "estiva/mesh.h"
#include "estiva/ary.h"
#include "estiva/op.h"
#include "estiva/mx.h"
#include "estiva/foreach.h"
#include "estiva/que.h"


void nsA(MX **Ap, xyc *Z, nde *N)
{
  static double *S, *U, *V;
  static MX *M, *AX, *AY, *D, *HX, *HY;
  static MX *A;

  double tau = 0.1;
  long i, j, m, n;
    
  ary1(U,10000);
  for ( i = 0; i< 9999; i++ ) U[i] = 1.0;
  ary1(V,10000);
  for ( i = 0; i< 9999; i++ ) V[i] = 1.0;

  femdelta(S,Z,N);
  nsM(M,S,N);
  nsAX(AX,U,S,Z,N);
  nsAY(AY,V,S,Z,N);
  nsD(D,S,Z,N);
  nsHX(HX,S,Z,N);
  nsHY(HY,S,Z,N);

  m = dimp2(N);
  n = dim1(Z);
  initmx(*Ap,2*m+n+1,50);
  A = *Ap;
  
  for(i=1; i<=m; i++) for(j=1; j<=n; j++) mx(A,i,2*m+j)   = -mx(HX,i,j);
  for(i=1; i<=m; i++) for(j=1; j<=n; j++) mx(A,m+i,2*m+j) = -mx(HY,i,j);
  for(i=1; i<=n; i++) for(j=1; j<=m; j++) mx(A,2*m+i,j)   =  mx(HX,j,i);
  for(i=1; i<=n; i++) for(j=1; j<=m; j++) mx(A,2*m+i,m+j) =  mx(HY,j,i);

  for(i=1; i<=m; i++) for(j=1; j<=m; j++) {
      mx(A,i  , j  ) = mx(M,i,j)/tau + mx(D,i,j);
    }
  for(i=1; i<=m; i++) for(j=1; j<=m; j++) {
      mx(A,i+m, j+m) = mx(M,i,j)/tau + mx(D,i,j);
    }
}


int main(int argc, char **argv)
{
  static xyc *Z; static nde *N;
  static MX *A;
  static double *x;
  long i;
  
  initop(argc,argv);
  rectmesh(Z,N);

  ary1(x,dimp2(N)*2+dim1(Z)+1);
  for (i=1;i<=dimp2(N);i++) x[i] = -0.05;
  for (i=dimp2(N)+1;i<=dimp2(N)*2;i++) x[i] = -0.05;

  nsA(&A,Z,N);
  
  forgammap2(i,"zero" ,Z,N) { x[i] = 0.0;  x[dimp2(N)+i] = 0.0; }
  forgammap2(i,"gamma",Z,N) { x[i] = 0.05; x[dimp2(N)+i] = 0.0; }
  
  pltp2(x,Z,N);
  sleep(60);
  return 0;
}
