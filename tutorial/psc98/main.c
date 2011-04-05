#include <stdio.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/op.h"
#include "estiva/solver.h"
#include "estiva/std.h"

static void fprintmx(FILE *fp,MX *A)
{
  long i, j;
  fornonzeromx(A,i,j) fprintf(fp,"%ld %ld %e\n",i,j,mx(A,i,j));
  fclose(fp);
}

int main(int argc, char **argv){
  static MX *A;
  static int *JA;
  static double *AA, B, *b, *x;
  int i, j, n, w;
  initop(argc,argv);
  ary1(JA,10);
  ary1(AA,10);

  genmat(-1,JA,AA,&B);

  n = JA[0]; w = JA[2];
  fprintf(stderr,"n=%d, JA[1]=%d, w=%d B=%e\n",n,JA[1],w,B);
  initmx(A,n+1,w+1);
  ary1(b,n+1);
  ary1(x,n+1);
  for ( i=1; i<=n; i++) {
    forall (0,j,w-1) {
      JA[j] =  -1;
      AA[j] = 0.0;
    }
    genmat(i,JA,AA,&B);
    for ( j=0; j<w; j++) if (JA[j] != -1){
	if( JA[j] <=0 || n < JA[j] ) {
	  ;
	} else {
	  mx(A,i,JA[j]) = AA[j]; 
	  b[i] = B;
	}
    }
  }
  fprintmx(fopen("foo","w"),A);

  if (!defop("--mpi"))
    mpisolver(A,x,b);
  else
    solver(A,x,b);

  chkval(stdout,n,&x[1]);
  sendcommand(999);
}
