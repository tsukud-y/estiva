#include <stdio.h>
#include <stdlib.h>
#include <estiva/op.h>
#include <estiva/ary.h>
#include <estiva/mesh.h>


static void estiva_p1plt(FILE *fp, xyc *Z, nde *N, double *u)
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

static void conta(FILE *fp,double *P, xyc *Z, nde *N, double k)
{
  long e, a, b, c;
  double ax,ay,az, bx,by,bz, cx,cy,cz, x0,y0, x1,y1;

  for(e=1;e<=dim1(N);e++){
    a = N[e].a, b = N[e].b, c = N[e].c;

    ax = Z[a].x, ay = Z[a].y, az = P[a];
    bx = Z[b].x, by = Z[b].y, bz = P[b];
    cx = Z[c].x, cy = Z[c].y, cz = P[c];


    if((az>k&&bz<k&&cz<k)||(az<k&&bz>k&&cz>k)){
      x0 = (k-az)*(bx-ax)/(bz-az) + ax;
      y0 = (k-az)*(by-ay)/(bz-az) + ay;
      x1 = (k-az)*(cx-ax)/(cz-az) + ax;
      y1 = (k-az)*(cy-ay)/(cz-az) + ay;
      fprintf(fp,"%f %f %f \n",x0,y0,k);
      fprintf(fp,"%f %f %f \n",x1,y1,k);
      fprintf(fp,"\n\n");
    }
    if((bz>k&&cz<k&&az<k)||(bz<k&&cz>k&&az>k)){
      x0 = (k-bz)*(cx-bx)/(cz-bz) + bx;
      y0 = (k-bz)*(cy-by)/(cz-bz) + by;
      x1 = (k-bz)*(ax-bx)/(az-bz) + bx;
      y1 = (k-bz)*(ay-by)/(az-bz) + by;
      fprintf(fp,"%f %f %f \n",x0,y0,k);
      fprintf(fp,"%f %f %f \n",x1,y1,k);
      fprintf(fp,"\n\n");
    }
    if((cz>k&&az<k&&bz<k)||(cz<k&&az>k&&bz>k)){
      x0 = (k-cz)*(ax-cx)/(az-cz)+cx;
      y0 = (k-cz)*(ay-cy)/(az-cz)+cy;
      x1 = (k-cz)*(bx-cx)/(bz-cz)+cx;
      y1 = (k-cz)*(by-cy)/(bz-cz)+cy;
      fprintf(fp,"%f %f %f\n",x0,y0,k);
      fprintf(fp,"%f %f %f\n",x1,y1,k);
      fprintf(fp,"\n\n");
    }
  }
}

static void estiva_contaplt(FILE *fp, xyc *Z, nde *N, double *u)
{
  if ( fp == NULL ) {
    fp = fopen("/tmp/plt.tmp","w");
    {
      double x, h, umax = u[1], umin = u[1];
      int i;

      for (i=2; i<=dim1(Z); i++) {
	if(umax < u[i]) umax = u[i];
	if(umin > u[i]) umin = u[i];
      }
      h = (umax-umin)/20.0;

      if( 0.0 != atof(getop("-conta")) ) 
	h = (umax-umin)/(1.0+atof(getop("-conta")));
      
      for (x = umin; x<=umax; x+=h) conta(fp,u,Z,N,x);

      fclose(fp);
      fp = popen("gnuplot","w");
      fprintf(fp,"splot '/tmp/plt.tmp' w l\n");
      fflush(fp);
    }
  }
}


void estiva_plt(FILE *fp, xyc *Z, nde *N, double *u)
{
  if (defop("-conta")) {
    estiva_contaplt(fp, Z, N, u);
    return ;
  }
  estiva_p1plt(fp, Z, N, u);
}

