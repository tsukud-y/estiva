#include <stdio.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/solver.h>

main(){
  int i, n = 10000;
  static MX *A; static double *x, *b;
  initmx(A,n+1,8); ary1(x, n+1); ary1(b, n+1);

  for ( i =1; i<=n; i++){
    mx(A,i,i) = 1.0; b[i] = 1.0;
  }
  solver(A,x,b);

  printf("%f\n", x[1] );
  printf("%f\n", x[2] );
}
