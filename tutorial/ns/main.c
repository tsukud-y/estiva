#include "ns.h"
#include "fem.h"


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
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
    b[i]          = 100.0;
    b[m+i]        = 0.0; 
  }

  i = NUM;
  for(j=1; j<=NUM; j++) mx(A,i,j) = 0.0;
  mx(A,i,i) = 1.0; 
  b[i] = 1.0; 
}


void b_(double *b, xyc *Z, nde *N, MX *M, double *x, double t)
{
  static double *Ux, *Uy, *bx, *by;
  long i, m, n;

  m = dimp2(N); 
  n = dim1(Z); 

  ary1(Ux,m+1);
  ary1(Uy,m+1);

  for ( i = 1; i <= m; i++ ) {
    Ux[i] = x[i];
    Uy[i] = x[i+m];
  }

  mulmx(bx,M,Ux);
  mulmx(by,M,Uy);

  for ( i = 1; i <= m; i++ ) {
    b[i  ] = t*bx[i];
    b[i+m] = t*by[i];
  }

  for ( i = 1; i <= n; i++ ) b[i + 2*m] = 0.0;

}


void nsA(MX **Ap, double *x, double *b, xyc *Z, nde *N, MX *K, MX *M, MX *Hx, MX *Hy, MX *AX, MX *AY, double t)
{
  static MX *A;
  long i, j, NUM, m, n;

  m   = dimp2(N); 
  n   = dim1(Z); 
  NUM = m*2+n;
  initmx(*Ap, NUM+1, 50); A = *Ap; 

  for ( i = 1; i <= m; i++ ) for ( j = 1; j <= m; j++ ) {
      mx(A,  i,   j) = mx(M,i,j) + t*mx(K,i,j) + t*mx(AX,i,j);
      mx(A,m+i, m+j) = mx(M,i,j) + t*mx(K,i,j) + t*mx(AY,i,j);
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
  static MX *A, *K, *M, *Hx, *Hy, *AX, *AY; static double *x, *b, *S;
  long  i, k, kn = 1, m, n, NUM;
  double t = 0.001;

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

  
  if ( defop("-kn") ) kn = atoi(getop("-kn"));
  for ( k = 1; k <= kn; k++ ) {
    nsAX(AX,x,S,Z,N);
    nsAY(AY,x+m,S,Z,N);
    nsA(&A,x,b,Z,N,K,M,Hx,Hy,AX,AY,t);
    b_(b,Z,N,M,x,t);
    boundary_condition(Z,N,A,b);
    solver(A,x,b);
   }
  pltp2(x,Z,N);
  sleep(30);
  return 0;
}
