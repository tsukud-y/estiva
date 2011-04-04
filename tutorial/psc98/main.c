#include <stdio.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/op.h"
#include "estiva/solver.h"
#include "estiva/std.h"

int main(int argc, char **argv){
  static MX *A;
  static long *JA;
  static double *AA, B, *b, *x;
  int i, j, n;
  initop(argc,argv);
  ary1(JA,9);
  ary1(AA,9);

  genmat(-1,JA,AA,&B);

  n = JA[0];
  initmx(A,n+1,9);
  ary1(b,n+1);
  ary1(x,n+1);
  for ( i=1; i<=n; i++) {
    genmat(i,JA,AA,&B);
    for ( j=0; j<9; j++) if (JA[j] != -1){
	if( JA[j] < -1 || n < JA[j] ) {
	} else {
	  printf("%ld\n",JA[j]);
	  mx(A,i,JA[j]) = AA[j]; 
	  b[i] = B;
	}
    }
  }
  printf("hello\n");

  if (!defop("--mpi"))
    mpisolver(A,x,b);
  else
    solver(A,x,b);

  chkval(stdout,n,&x[1]);
  sendcommand(999);
}
