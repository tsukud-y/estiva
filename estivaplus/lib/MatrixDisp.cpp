#include "estivaplus.h"
#include <unistd.h>

void MatrixDisp(Matrix &A)
{
  FILE *pp;
  pp = popen("gnuplot","w");
  fprintf(pp,"set yrange [%ld:0]\n",A.capacity());
  forMatrixNonzero(A) {
    fprintf(pp,"set label \"%.1f\" at %u, %u;\n",A[i][j],j+1,i+1);
  }
  fprintf(pp,"plot '-' title \"\"\n");
  forMatrixNonzero(A) {
    fprintf(pp,"%u %u\n",j+1,i+1);
  }
  fprintf(pp,"e\n");
  fflush(pp);
  sleep(60*60);
  exit(0);
}
