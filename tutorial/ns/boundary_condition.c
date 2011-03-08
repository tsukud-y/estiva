#include "ns.h"


void estiva_boundary_condition(MX *A, double *b)
{
  static double velocity = 0.1;
  long i, j, NUM, m, n;
  static xyc *Z; static nde *N; static double *S;

  getZNS(Z,N,S);
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
  b[i] = 0.001;
}
