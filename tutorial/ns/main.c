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
#include "estiva/solver.h"


void nsA(MX **Ap, double *x, double *b, xyc *Z, nde *N)
{
  static double *S;
  static MX *M, *K, *Hx, *Hy;
  double t = 0.001;

  femdelta(S,Z,N);
  nsM(M,S,N);
  nsD(K,S,Z,N);
  nsHX(Hx,S,Z,N);
  nsHY(Hy,S,Z,N);


  MX *A; long i, j, NUM, m, n;
  m = dimp2(N);
  n = dim1(Z);

  NUM = m*2+n;
  initmx(*Ap, NUM+1, 50);  
  A = *Ap;

  for (i=1;i<=NUM;i++) for (j=1;j<=NUM;j++) mx(A,i,j) = 0.0;


  for(i=1;i<=m;i++) for(j=1; j<=m; j++){
      mx(A,  i,   j) = mx(M,i,j) + t*mx(K,i,j);
      mx(A,m+i, m+j) = mx(M,i,j) + t*mx(K,i,j);
    }

  for(i=1;i<=m;i++) for(j=1; j<=n; j++){
      mx(A,i,2*m+j) = -t*mx(Hx,i,j);
      mx(A,2*m+j,i) = -t*mx(Hx,i,j);
    }

  for(i=1;i<=m;i++) for(j=1; j<=n; j++){
      mx(A,m+i,2*m+j) = -t*mx(Hy,i,j);
      mx(A,2*m+j,m+i) = -t*mx(Hy,i,j);
    }

  forgammap2(i,"zero",Z,N) { 
    for (j=1; j<=NUM; j++) {
      mx(A,i,j)   = 0.0; 
      mx(A,m+i,j) = 0.0;
    }
    mx(A,i,i)     = 1.0; 
    mx(A,m+i,m+i) = 1.0;
    b[i]   = 0.0;  
    b[m+i] = 0.0; 
  }
  forgammap2(i,"gamma",Z,N) { 
    for (j=1; j<=NUM; j++) {
      mx(A,i,j)   = 0.0; 
      mx(A,m+i,j) = 0.0;
    }
    mx(A,i,i)     = 1.0; 
    mx(A,m+i,m+i) = 1.0;
    b[i]        = 0.0625;
    b[m+i]      = 0.0; 
  }

  for(j=1; j<=NUM; j++) mx(A,NUM,j) = 0.0;
  mx(A,NUM,NUM) = 1.0; 
  b[NUM] = 1.0; 
}


int main(int argc, char **argv)
{
  static xyc *Z; static nde *N; static MX *A; static double *x, *b;
  long  m, n, NUM;
  
  initop(argc,argv);
  rectmesh(Z,N);

  m = dimp2(N);
  n = dim1(Z);
  NUM = m*2+n;
  ary1(x,NUM+1); ary1(b,NUM+1);
  nsA(&A,x,b,Z,N);
  solver(A,x,b);
  pltp2(x,Z,N);
  sleep(60);
  return 0;
}
