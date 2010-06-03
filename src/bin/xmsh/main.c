#include <stdio.h>
#include <estiva/ary.h>
#include <estiva/op.h>
#include <estiva/mesh.h>


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


int main(int argc, char **argv){
  static xyc *Z;
  static nde *N;

  initop(argc, argv);

  fp2mesh(stdfp(),Z, N);

  pltmsh(stdout,Z,N);
  return 0;
}
