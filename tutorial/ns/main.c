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


double  mij(long i, long j);
double  kij(long i, long j);
double hxij(long i, long j);
double hyij(long i, long j);


void boundary_condition(xyc *Z, nde *N, MX *A, double *b)
{
  long i, j, NUM, m, n;
  m = dimp2(N); n = dim1(Z); NUM = m*2+n;

  for ( i = 1; i <= n; i++ ) if ( Z[i].label && !strcmp(Z[i].label,"zero") ) {
      for(j=1; j<=NUM; j++) { mx(A,i,j) = 0.0; mx(A,m+i,j) = 0.0; } 
      mx(A,i,i)     = 1.0;
      mx(A,m+i,m+i) = 1.0;
      b[i]          = 0.0;  
      b[i+m]        = 0.0;
    }


  forgammap2(i,"south",Z,N) {
    for(j=1; j<=NUM; j++) { mx(A,i+m,j) = 0.0; }
    mx(A,i+m,i+m) = 1.0;
    b[i+m]        = 0.0;  
  }

  forgammap2(i,"east",Z,N) {
    for(j=1; j<=NUM; j++) { mx(A,i,j) = 0.0;} 
    mx(A,i,i) = 1.0;
    b[i]      = 0.0; 
  }

  forgammap2(i,"west",Z,N) {
    for(j=1; j<=NUM; j++) { mx(A,i,j) = 0.0; }
    mx(A,i,i) = 1.0;
    b[i]      = 0.0; 
  }

  forgammap2(i,"north",Z,N) {
    for(j=1; j<=NUM; j++) { mx(A,i,j) = 0.0; mx(A,m+i,j) = 0.0; }
    mx(A,i,i)     = 1.0;
    mx(A,m+i,m+i) = 1.0;
    b[i]          = 0.1;
    b[m+i]        = 0.0; 
  }

  i = NUM;
  for(j=1; j<=NUM; j++) mx(A,i,j) = 0.0;
  mx(A,i,i) = 1.0; 
  b[i] = 1.0; 
}


void b_(double *b, xyc *Z, nde *N, MX *M, double *x)
{
  static double *bx, *xx, *by, *xy;
  long i, j, m, n;

  m = dimp2(N); n = dim1(Z); 

  ary1(xx,m+1);
  ary1(xy,m+1);

  
  for ( i = 1; i <= 2*m+n; i++ ) b[i] = 0.0;

  for ( i = 1; i <= m; i++ ) {
    xx[i] = x[i];
    xy[i] = x[i+m];
  }
  mulmx(bx,M,xx);
  mulmx(by,M,xy);
  for ( i = 1; i <= m; i++ ) {
    b[i]   = bx[i];
    b[i+m] = by[i];
  }
}


void nsA(MX **Ap, double *x, double *b, xyc *Z, nde *N, MX *K, MX *M, MX *Hx, MX *Hy)
{
  static MX *A;
  long i, j, NUM, m, n;
  double t = 0.001, Re = 1.0;

  m   = dimp2(N); 
  n   = dim1(Z); 
  NUM = m*2+n;
  initmx(*Ap, NUM+1, 50); A = *Ap; 

  for ( i = 1; i <= m; i++ ) for ( j = 1; j <= m; j++ ) {
      mx(A,  i,   j) = mx(M,i,j) + t*mx(K,i,j);
      mx(A,m+i, m+j) = mx(M,i,j) + t*mx(K,i,j);
    }
  for ( i = 1; i <= m; i++ ) for ( j = 1; j <= n; j++ ) {
      mx(A,    i,2*m+j) = -t*mx(Hx,i,j);
      mx(A,2*m+j,    i) = -t*mx(Hx,i,j);
      mx(A,  m+i,2*m+j) = -t*mx(Hy,i,j);
      mx(A,2*m+j,  m+i) = -t*mx(Hy,i,j);
    }
}


int main(int argc, char **argv)
{
  static xyc *Z; static nde *N; 
  static MX *A, *K, *M, *Hx, *Hy; static double *x, *b, *S;
  long  i, k, m, n, NUM;
  
  initop(argc,argv);
  rectmesh(Z,N);

  m   = dimp2(N); 
  n   = dim1(Z); 
  NUM = m*2 + n;

  ary1(x,NUM+1); 
  ary1(b,NUM+1);

  femdelta(S,Z,N); 
  setZNS(Z,N,S);
  genP2P2mx(&M,mij);
  genP2P2mx(&K,kij);
  genP2P1mx(&Hx,hxij);
  genP2P1mx(&Hy,hyij);
  nsA(&A,x,b,Z,N,K,M,Hx,Hy);

  for ( k = 1; k <= 10; k++ ) {
    b_(b,Z,N,M,x);
    boundary_condition(Z,N,A,b);
    solver(A,x,b);
    for ( i = 1; i <= m * 2; i++ ) x[i] *= 40.0;
    pltp2(x,Z,N);
    for ( i = 1; i <= m * 2; i++ ) x[i] /= 40.0;
   }
  sleep(30);
  return 0;
}
