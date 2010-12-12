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
  static double velocity = 0.1;
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

  if ( defop("-velocity") ) velocity = atof(getop("-velocity"));

  forgammap2(i,"north",Z,N) {
    for(j=1; j<=NUM; j++) { mx(A,i,j) = 0.0; mx(A,m+i,j) = 0.0; }
    mx(A,i,i)     = 1.0;
    mx(A,m+i,m+i) = 1.0;
    b[i]          = velocity;
    b[m+i]        = 0.0; 
  }

  i = NUM;
  for(j=1; j<=NUM; j++) mx(A,i,j) = 0.0;
  mx(A,i,i) = 1.0; 
  b[i] = 1.0; 
}


void nsRhs(double *b, xyc *Z, nde *N, MX *M, double *x, double t)
{
  long   i, j, NUM, m, n;

  m = dimp2(N);
  n = dim1(Z);

  NUM = m*2+n;

  for( i = 1; i <= NUM; i++ ) b[i] = 0.0;

  for( i = 1; i <= m; i++ ) for ( j = 1; j <= m; j++ ) {
      b[  i] += mx(M,i,j)*x[  j];
      b[m+i] += mx(M,i,j)*x[m+j];
    }


  /*
  static int init = 0;
  static MX *M2;

  if ( !init ) {
    long i, j, m, n;
    init = 1;
    m = dimp2(N); 
    n = dim1(Z); 

    initmx(M2, 2*m+n, 50);
    
    for ( i = 1; i <= m; i++ ) for ( j = 1; j <=m; j++ ) {
	mx(M2,  i,  j) = mx(M,i,j);
	mx(M2,m+i,m+j) = mx(M,i,j);
      }
  }
  mulmx(b,M2,x);
  */
}


void nsA(MX **Ap, double *x, double *b, xyc *Z, nde *N, MX *K, MX *M, MX *Hx, MX *Hy, MX *AX, MX *AY, double tau)
{
  static MX *A;
  long   i, j, NUM, m, n;

  m   = dimp2(N); 
  n   = dim1(Z); 
  NUM = m*2+n;
  initmx(*Ap, NUM+1, 50); A = *Ap; 

  for ( i = 1; i <= m; i++ ) for ( j = 1; j <= m; j++ ) {
      mx(A,  i,   j) = mx(M,i,j) + tau*mx(K,i,j) + tau*mx(AX,i,j);
      mx(A,m+i, m+j) = mx(M,i,j) + tau*mx(K,i,j) + tau*mx(AY,i,j);
    }
  for ( i = 1; i <= m; i++ ) for ( j = 1; j <= n; j++ ) {
      mx(A,    i,2*m+j) = -tau*mx(Hx,i,j);
      mx(A,2*m+j,    i) = -tau*mx(Hx,i,j);
      mx(A,  m+i,2*m+j) = -tau*mx(Hy,i,j);
      mx(A,2*m+j,  m+i) = -tau*mx(Hy,i,j);
    }
}


int main(int argc, char **argv)
{
  static xyc *Z; static nde *N; 
  static MX *A, *K, *M, *Hx, *Hy, *AX, *AY; static double *x, *b, *S;
  long   k, kn = 1, m, n, NUM;
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
    nsRhs(b,Z,N,M,x,t);
    boundary_condition(Z,N,A,b);
    solver(A,x,b);
    printf("kn - k = %ld\n",kn-k);
    pltp2(x,Z,N);
   }
  return 0;
}
