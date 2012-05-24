#include "estivaplus.h"
#include "stwart.h"

void printMatrix(Matrix &B2)
{
  long i, j, M=B2.size();
  for (i=0;i<M;i++){
    for (j=0;j<M;j++) 
      if (B2[i][j] == 0.0) printf("  ");
      else printf("%ld ",(long)B2[i][j]);
    printf("\n");
  }
}

