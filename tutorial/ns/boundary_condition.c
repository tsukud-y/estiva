#include "ns.h"


void estiva_boundary_condition(MX *A, double *b)
{
  static double velocity = 1000.0;
  long i, NUM, m, n;
  static xyc *Z; static nde *N; static double *S;

  getZNS(Z,N,S);
  m = dimp2(N); n = dim1(Z); NUM = m*2+n;

  for ( i = 1; i <= n; i++ ) if ( Z[i].label && !strcmp(Z[i].label,"zero") ) {
      zerofillrow(A,i);
      zerofillrow(A,m+i);
      mx(A,i,i)     = 1.0;
      mx(A,m+i,m+i) = 1.0;
      b[i]          = 0.0;
      b[i+m]        = 0.0;
    }


  forgammap2(i,"south",Z,N) {
    zerofillrow(A,i);
    zerofillrow(A,i+m);
    mx(A,i,i)     = 1.0;
    mx(A,i+m,i+m) = 1.0;
    b[i]          = 0.0;
    b[i+m]        = 0.0;
  }

  forgammap2(i,"east",Z,N) {
    zerofillrow(A,i);
    zerofillrow(A,i+m);
    mx(A,i,i)     = 1.0;
    mx(A,i+m,i+m) = 1.0;
    b[i]          = 0.0;
    b[i+m]        = 0.0;
  }

  forgammap2(i,"west",Z,N) {
    zerofillrow(A,i);
    zerofillrow(A,i+m);
    mx(A,i,i)     = 1.0;
    mx(A,i+m,i+m) = 1.0;
    b[i]          = 0.0;
    b[i+m]        = 0.0;
  }

  if ( defop("-velocity") ) velocity = atof(getop("-velocity"));

  forgammap2(i,"north",Z,N) {
    zerofillrow(A,i);
    zerofillrow(A,m+i);
    mx(A,i,i)     = 1.0;
    mx(A,m+i,m+i) = 1.0;
    b[i]          = velocity;
    b[m+i]        = 0.0;
  }

  i = NUM;
  zerofillrow(A,i);
  mx(A,i,i) = 1.0;
  b[i] = 0.001;
}
