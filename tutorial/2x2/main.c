#include <stdio.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/solver.h>

main(){
  static MX *A; static double *x, *b;
  initmx(A,3,3); ary1(x, 3); ary1(b, 3);

  mx(A,1,1) = 1.0; mx(A,1,2) = 0.0; b[1] = 3.0;
  mx(A,2,1) = 0.0; mx(A,2,2) = 1.0; b[2] = 2.0;
  solver(A,x,b);

  printf("%f\n", x[1] );
  printf("%f\n", x[2] );
}
