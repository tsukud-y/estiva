#include <stdio.h>
#include <estiva/op.h>
#include <estiva/ary.h>
#include <estiva/mesh.h>


void estiva_plt(FILE *fp, xyc *Z, nde *N, double *u)
{
  if ( fp == NULL ) {
    fp = fopen("/tmp/plt.tmp","w");
    {
      long e, a, b, c;
      for(e=1;e<=dim1(N);e++){
	a = N[e].a, b = N[e].b, c = N[e].c;
	fprintf(fp,"%f %f %f\n",Z[a].x,Z[a].y,u[a]);
	fprintf(fp,"%f %f %f\n",Z[b].x,Z[b].y,u[b]);
	fprintf(fp,"%f %f %f\n",Z[c].x,Z[c].y,u[c]);
	fprintf(fp,"%f %f %f\n",Z[a].x,Z[a].y,u[a]);
	fprintf(fp,"\n\n");
      }
      fclose(fp);
      fp = popen("gnuplot","w");
      fprintf(fp,"splot '/tmp/plt.tmp' w l\n");
      fflush(fp);
    }
  }
}
