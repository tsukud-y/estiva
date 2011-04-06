#include <stdio.h>
#include <stdlib.h>
#include "estiva/mx.h"
#include <sys/types.h>
#include <unistd.h>


void estiva_pltmx(MX *A)
{
  if ( fork() == 0 ) {
    FILE *pp;
    long i, j;

    pp = popen("gnuplot","w");
    fprintf(pp,"set yrange [%ld:0]\n",A->m+1);

    mx(A,1,1) = mx(A,1,1);
    for(i=0; i<A->m; i++)
      for (j=0; j<A->w; j++) 
	if(A->A[i][j] != 0.0) { 
	  fprintf(pp,"set label \"%.1f\" at %ld, %ld;\n",A->A[i][j],A->IA[i][j],i+1);
	}

    fprintf(pp,"plot '-' title \"\"\n");
    mx(A,1,1) = mx(A,1,1);
    for(i=0; i<A->m; i++)
      for (j=0; j<A->w; j++) 
	if(A->A[i][j] != 0.0) { 
	  fprintf(pp,"%ld %ld\n",A->IA[i][j],i+1);
	}
    fprintf(pp,"e\n");
    fflush(pp);
    sleep(60*3);
    pclose(pp);
    exit(0);
  }
}
