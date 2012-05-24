#include "estivaplus.h"

void BoundaryCondition(Matrix &A, Vector &b)
{
  static xyc *Z; static nde *N; static double *S;
  long i, m;

  estiva_getZNS(&Z,&N,&S);
  m = dimp2(N);

  forgammap1(i,"zero",Z) {
    A[i-1].clear();
    A[m+i-1].clear();
    A[i-1][i-1]     = 1.0;
    A[m+i-1][m+i-1] = 1.0;
    b[i-1]          = 0.0;
    b[i+m-1]        = 0.0;
  }
  
  forgammap2(i,"cylinder",Z,N) {
    A[i-1].clear();
    A[m+i-1].clear();
    A[i-1][i-1]     = 1.0;
    A[m+i-1][m+i-1] = 1.0;
    b[i-1]          = 0.0;
    b[i+m-1]        = 0.0;
  }
  
  forgammap2(i,"south",Z,N) {
    A[i+m-1].clear();
    A[i+m-1][i+m-1] = 1.0;
    b[i+m-1]        = 0.0;
  }
  
  forgammap2(i,"west",Z,N) {
    A[i-1].clear();
    A[i+m-1].clear();
    A[i-1][i-1]     = 1.0;
    A[i+m-1][i+m-1] = 1.0;
    b[i-1]          = 1.0;
    b[i+m-1]        = 0.0;
  }

  forgammap2(i,"north",Z,N) {
    A[m+i-1].clear();
    A[m+i-1][m+i-1] = 1.0;
    b[m+i-1]        = 0.0;
  }


#if 0
  long j;
  for(i=0;i<200;i++) {
    for (j=2*m; j<2*m+200;j++)
      if(A[i][j] == 0) printf(" ");
      else             printf("%f",A[i][j]);
    printf("|\n");
  }
#endif

#if 0  
  long NUM = b.size();


  for(i=0;i<NUM;i++) printf("%ld,%f\n",i, A[i][i]);
  printf("m = %ld\n",m);


  long NUM = b.size();
  A[NUM-1].clear();
  A[NUM-1][NUM-1] = 1.0;
  b[NUM-1] = 0.0;
#endif

}
