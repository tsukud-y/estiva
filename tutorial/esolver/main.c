#include <stdio.h>
#include <estiva/ary.h>
#include <estiva/esolver.h>

main(){
  static double **A, *x;

  ary2(A,3+1,3+1); ary1(x,3+1);

  A[1][1] = 1.0;  A[1][2] = 0.0;  A[1][3] = 0.0;
  A[2][1] = 0.0;  A[2][2] = 1.0;  A[2][3] = 0.0;
  A[3][1] = 0.0;  A[3][2] = 0.0;  A[3][3] = 1.0;

  x[1] = 1.0;
  x[2] = 1.0;  
  x[3] = 1.0;

  printf("eigenvalue = %f\n",esolver(A,x));
  printf("%f %f %f\n",x[1],x[2],x[3]);
}
