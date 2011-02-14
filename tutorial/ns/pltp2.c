#include "fem.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva/mesh.h"
#include "estiva/foreach.h"
#include "estiva/ary.h"
#include "estiva/op.h"

static char *num2str(long n){
  long i;
  static char ret[100];
  sprintf(ret,"00000000%ld",n);
  i = strlen( ret );
  
  return &ret[i-8];
}

static void pltmsh(FILE *fp, xyc *Z, nde *N)
{

  long e, a, b, c, dim1N;
  dim1N = dim1(N);
  for(e=1;e<=dim1N;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    fprintf(fp,"%f %f\n",Z[a].x,Z[a].y);
    fprintf(fp,"%f %f\n",Z[b].x,Z[b].y);
    fprintf(fp,"%f %f\n",Z[c].x,Z[c].y);
    fprintf(fp,"%f %f\n",Z[a].x,Z[a].y);
    fprintf(fp,"\n\n");
  }
}

void estiva_pltp2(double *x, xyc * Z, nde *N)
{
  FILE *pp;
  static double *u, *v;
  long e, m, dim1N;
  int n;
  double x0, y0;
  double scale = 1.0;
  static long sequence = 0;
  static char filename[400];
  long arrow = 1;

  sprintf(filename,"anime/%s.gnuplot",num2str(sequence++));
  pp = fopen(filename,"w");


  u = x;
  v = x + dimp2(N);

  if ( defop("-plotscale") ) scale = atof(getop("-plotscale"));

  dim1N = dim1(N);
  for (e=1; e<=dim1N; e++)
    foreach(n) &N[e].a, &N[e].b, &N[e].c, end {
      fprintf(pp,"set arrow %ld from %f,%f to %f,%f\n",arrow++,Z[n].x,Z[n].y,Z[n].x+u[n]*scale,Z[n].y+v[n]*scale);
      if ( Z[n].label != NULL ) 
	fprintf(pp,"set label \"%s\" at %f, %f;\n",Z[n].label, Z[n].x, Z[n].y);
    }



  for (e=1; e<=dim1N; e++) {
    static char *label="zero";
    m = N[e].A;
    x0 = (Z[N[e].b].x + Z[N[e].c].x)/2.0, y0 = (Z[N[e].b].y + Z[N[e].c].y)/2.0;
    fprintf(pp,"set arrow %ld from %f,%f to %f,%f\n",arrow++,x0,y0,x0+u[m]*scale,y0+v[m]*scale);
    if ( (Z[N[e].b].label && !strcmp(Z[N[e].b].label, label)) ||

         (Z[N[e].c].label && !strcmp(Z[N[e].c].label, label))  ) {
      if (Z[N[e].b].label && Z[N[e].c].label) 
	fprintf(pp,"set label \"%s\" at %f, %f;\n",label,x0,y0);
    }
    m = N[e].B;
    x0 = (Z[N[e].c].x + Z[N[e].a].x)/2.0, y0 = (Z[N[e].c].y + Z[N[e].a].y)/2.0;
    fprintf(pp,"set arrow %ld from %f,%f to %f,%f\n",arrow++,x0,y0,x0+u[m]*scale,y0+v[m]*scale);
    if ( (Z[N[e].c].label && !strcmp(Z[N[e].c].label, label)) ||
         (Z[N[e].a].label && !strcmp(Z[N[e].a].label, label))  ) {
      if (Z[N[e].c].label && Z[N[e].a].label)
	fprintf(pp,"set label \"%s\" at %f, %f;\n",label,x0,y0);
    }

    m = N[e].C;
    x0 = (Z[N[e].a].x + Z[N[e].b].x)/2.0, y0 = (Z[N[e].a].y + Z[N[e].b].y)/2.0;
    fprintf(pp,"set arrow %ld from %f,%f to %f,%f\n",arrow++,x0,y0,x0+u[m]*scale,y0+v[m]*scale);
    if ( (Z[N[e].a].label && !strcmp(Z[N[e].a].label, label)) ||
         (Z[N[e].b].label && !strcmp(Z[N[e].b].label, label))  ) {
      if (Z[N[e].a].label && Z[N[e].b].label)
	fprintf(pp,"set label \"%s\" at %f, %f;\n",label,x0,y0);
    }
  }

  fprintf(pp,"plot '-' title \"\" with lines\n");
  pltmsh(pp,Z,N);
  fprintf(pp,"e\n");
  fflush(pp);
  fclose(pp);
}
