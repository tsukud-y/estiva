#include <stdio.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/solver.h>
#include <estiva/op.h>

main(int argc, char **argv){
  int i, n = 512;
  static MX *A; static double *x, *b;
  initop(argc,argv);
  initmx(A,n+1,8); ary1(x, n+1); ary1(b, n+1);

  for ( i =1; i<=n; i++){
    mx(A,i,i) = 2.0; b[i] = 1.0;
  }
  for (i=1;i<n;i++) mx(A,i,i+1) = -1.0;
  for (i=1;i<n;i++) mx(A,i+1,i) = -1.0;

  solver(A,x,b);

  for (i=1;i<=n;i++) printf("%f\n", x[i] );
}
