#include "ns.h"
#include "fem.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "estiva/mesh.h"
#include "estiva/ary.h"
#include "estiva/op.h"
#include "estiva/mx.h"
#include "estiva/foreach.h"
#include "estiva/que.h"
#include "estiva/solver.h"

extern double mij(long i, long j);
extern double kij(long i, long j);
extern double hxij(long i, long j);
extern double hyij(long i, long j);

void nsA(MX **Ap, double *x, double *b, xyc *Z, nde *N)
{
  static MX *A, *M, *K, *AX, *AY, *Hx, *Hy;
  static double *S, *U, *V, *u, *v;
  long i, j, NUM, m, n;
  double t = 0.001, Re = 1.0;

  m = dimp2(N); n = dim1(Z); NUM = m*2+n;

  initmx(*Ap, NUM+1, 50); A = *Ap; 
  ary1(U,m+1); ary1(V,m+1);

  femdelta(S,Z,N);
  setZNS(Z,N,S);
  genP2P2mx(&K,kij);
  genP2P2mx(&M,mij);
  genP2P1mx(&Hx,hxij);
  genP2P1mx(&Hy,hyij);
  for(i=1;i<=m;i++) U[i] = b[i];
  for(i=1;i<=m;i++) V[i] = b[i+m];
  nsAX(AX, U, S, Z, N);
  nsAY(AY, V, S, Z, N);
  mulmx(u,M,U);
  mulmx(v,M,V);

  for (i=1;i<=m;i++) b[i] = u[i];
  for (i=1;i<=m;i++) b[i+m] = v[i];

  for(i=1;i<=m;i++) for(j=1; j<=m; j++){
      mx(A,  i,   j) = mx(M,i,j)/t  + mx(K,i,j)/Re;
      mx(A,m+i, m+j) = mx(M,i,j)/t  + mx(K,i,j)/Re;
    }
  for(i=1;i<=m;i++) for(j=1; j<=n; j++){
      mx(A,i,2*m+j)   = -mx(Hx,i,j);
      mx(A,2*m+j,i)   = -mx(Hx,i,j);
      mx(A,m+i,2*m+j) = -mx(Hy,i,j);
      mx(A,2*m+j,m+i) = -mx(Hy,i,j);
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
    b[i]        = 1.0;
    b[m+i]      = 0.0; 
  }

  for(j=1; j<=NUM; j++) mx(A,NUM,j) = 0.0;
  mx(A,NUM,NUM) = 1.0; 
  b[NUM] = 1.0; 

}


int main(int argc, char **argv)
{
  static xyc *Z; static nde *N; static MX *A; static double *x, *b;
  long  i, k, m, n, NUM;
  
  initop(argc,argv);
  rectmesh(Z,N);
  m = dimp2(N);
  n = dim1(Z);
  NUM = m*2+n;

  ary1(x,NUM+1); ary1(b,NUM+1);

  for(k=1;k<=1;k++) {
    nsA(&A,x,b,Z,N);
    solver(A,x,b);
    pltp2(x,Z,N);
    for (i=1;i<=NUM;i++) b[i]=x[i];
  }
  sleep(60);
  return 0;
}
