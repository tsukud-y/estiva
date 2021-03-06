#include "ns.h"

void estiva_boundary_condition(MX *A, double *b)
{
  static double velocity = 1.0;
  long i, NUM, m, n;
  static xyc *Z; static nde *N; static double *S;
  
  getZNS(Z,N,S);
  m = dimp2(N); n = dim1(Z); NUM = m*2+n;

  forgammap1(i,"zero",Z) {
      zerofillrow(A,i);
      zerofillrow(A,m+i);
      mx(A,i,i)     = 1.0;
      mx(A,m+i,m+i) = 1.0;
      b[i]          = 0.0;
      b[i+m]        = 0.0;
    }


  forgammap2(i,"cylinder",Z,N) {
      zerofillrow(A,i);
      zerofillrow(A,m+i);
      mx(A,i,i)     = 1.0;
      mx(A,m+i,m+i) = 1.0;
      b[i]          = 0.0;
      b[i+m]        = 0.0;
    }


  forgammap2(i,"south",Z,N) {
    //zerofillrow(A,i);
    zerofillrow(A,i+m);
    //mx(A,i,i)     = 1.0;
    mx(A,i+m,i+m) = 1.0;
    //b[i]          = 0.0;
    b[i+m]        = 0.0;
  }

  forgammap2(i,"east",Z,N) {
    /*
    zerofillrow(A,i);
    zerofillrow(A,i+m);
    mx(A,i,i)     = 1.0;
    mx(A,i+m,i+m) = 1.0;
    b[i]          = 0.0;
    b[i+m]        = 0.0;
    */
  }

  forgammap2(i,"west",Z,N) {
    zerofillrow(A,i);
    zerofillrow(A,i+m);
    mx(A,i,i)     = 1.0;
    mx(A,i+m,i+m) = 1.0;
    b[i]          = 1.0;
    b[i+m]        = 0.0;
  }

  //if ( defop("-velocity") ) velocity = atof(getop("-velocity"));

  forgammap2(i,"north",Z,N) {
    //zerofillrow(A,i);
    zerofillrow(A,m+i);
    //mx(A,i,i)     = 1.0;
    mx(A,m+i,m+i) = 1.0;
    //b[i]          = 0.0;
    b[m+i]        = 0.0;
  }


#if 0
  //i = NUM-100;
  forgammap1(i,"zero",Z);
  i += 2*m;
  zerofillrow(A,i);
  mx(A,i,i) = 1.0;
  b[i] = 0.000;
#endif
}
