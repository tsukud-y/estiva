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


static void sleep_forever(void) 
{
  sleep(60*3);
}


void estiva_xmesh(xyc *Z)
{
  if ( fork() == 0 ) {
    long e, a, b, c, A, B, C;
    FILE *pp;
    static nde *N;
    delaunay(Z,N);
    pp = popen("gnuplot","w");

#if 0
    for (e=1; e<=dim1(N); e++) {
      a = N[e].a, b = N[e].b, c = N[e].c;
      A = N[e].A, B = N[e].B, C = N[e].C;
      fprintf(pp,"set label \"%ld\" at %f , %f\n",a,Z[a].x,Z[a].y);
      fprintf(pp,"set label \"%ld\" at %f , %f\n",b,Z[b].x,Z[c].y);
      fprintf(pp,"set label \"%ld\" at %f , %f\n",c,Z[b].x,Z[c].y);

      fprintf(pp,"set label \"%ld\" at %f , %f\n",A,(Z[b].x+Z[c].x)/2.0,(Z[b].y+Z[c].y)/2.0);
      fprintf(pp,"set label \"%ld\" at %f , %f\n",B,(Z[c].x+Z[a].x)/2.0,(Z[c].y+Z[a].y)/2.0);
      fprintf(pp,"set label \"%ld\" at %f , %f\n",C,(Z[a].x+Z[b].x)/2.0,(Z[a].y+Z[b].y)/2.0);

      fprintf(pp,"set label \"(%ld)\" at %f , %f\n",e,(Z[a].x+Z[b].x+Z[c].x)/3.0,(Z[a].y+Z[b].y+Z[c].y)/3.0);
    }
#endif
    for (e=1; e<=dim1(N); e++) {
      a = N[e].a, b = N[e].b, c = N[e].c;
      A = N[e].A, B = N[e].B, C = N[e].C;
      fprintf(pp,"set label \"%ld\" at %f , %f\n",a,Z[a].x,Z[a].y);
      fprintf(pp,"set label \"%ld\" at %f , %f\n",b,Z[b].x,Z[b].y);
      fprintf(pp,"set label \"%ld\" at %f , %f\n",c,Z[c].x,Z[c].y);
    }

    fprintf(pp,"plot '-' title \"\" with lines\n");
    pltmsh(pp,Z,N);
    fprintf(pp,"e\n");
    fflush(pp); 
    sleep_forever();
    pclose(pp);
    exit(0);
  }
}
