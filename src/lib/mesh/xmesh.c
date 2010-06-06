#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "estiva/op.h"
#include "estiva/mesh.h"
#include "estiva/ary.h"
#include <unistd.h>
#include <sys/types.h>


static void pltmsh(FILE *fp, xyc *Z, nde *N)
{
  long e, a, b, c;
  for(e=1;e<=dim1(N);e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    fprintf(fp,"%f %f\n",Z[a].x,Z[a].y);
    fprintf(fp,"%f %f\n",Z[b].x,Z[b].y);
    fprintf(fp,"%f %f\n",Z[c].x,Z[c].y);
    fprintf(fp,"%f %f\n",Z[a].x,Z[a].y);
    fprintf(fp,"\n\n");
  }
}


void sleep_forever(void) 
{
  sleep(60*3);
}


void estiva_xmesh(xyc *Z)
{
  if ( fork() == 0 ) {
    FILE *pp;
    static nde *N;
    delaunay(Z,N);
    pp = popen("gnuplot","w");
    fprintf(pp,"plot '-' title \"\" with lines\n");
    pltmsh(pp,Z,N);
    fprintf(pp,"e\n");
    fflush(pp); 
    sleep_forever();
    pclose(pp);
    exit(0);
  }
}
