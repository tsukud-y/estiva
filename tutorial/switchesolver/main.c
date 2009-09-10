#include <stdio.h>
#include <estiva/op.h>
#include <estiva/ary.h>
#include <estiva/esolver.h>

main(int argc, char **argv){
  static double **A, *x, lambda;     
  int i, n;

  initop(argc, argv);
  ary2(A,4,4); ary1(x,4);
  n = dim1(x);

  A[1][1] =  2.0;  A[1][2] = -1.0;  A[1][3] =  0.0;
  A[2][1] = -1.0;  A[2][2] =  2.0;  A[2][3] = -1.0;
  A[3][1] =  0.0;  A[3][2] = -1.0;  A[3][3] =  2.0;

  for (i=1; i<=n; i++ )
    x[i] = 1.0;

  lambda = esolver(A,x);

  printf("%s eigenvalue = %f\n",getop("-esolver"),lambda);
}

/*
実行(最大固有値を求める)
% ./a.out -esolver max

実行(最小固有値を求める)
% ./a.out -esolver min
*/
