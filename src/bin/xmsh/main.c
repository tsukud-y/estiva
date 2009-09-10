#include <stdio.h>
#include "msh.h"
#include "ary.h"
#include "op.h"

static void pltmsh(FILE *fp, xyc *Z, nde *N)
{
  long e, a, b, c;
  for(e=1;e<=dim1(N);e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    fprintf(fp,"%f %f %f\n",Z[a].x,Z[a].y);
    fprintf(fp,"%f %f %f\n",Z[b].x,Z[b].y);
    fprintf(fp,"%f %f %f\n",Z[c].x,Z[c].y);
    fprintf(fp,"%f %f %f\n",Z[a].x,Z[a].y);
    fprintf(fp,"\n\n");
  }
}


int main(int argc, char **argv){
  static xyc *Z;
  static nde *N;
  FILE *fp;
  int i;

  initop(argc, argv);
  fp = argf(argc, argv);

  fp2msh(fp,&Z, &N);

  pltmsh(stdout,Z,N);
}
