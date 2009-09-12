#include <stdio.h>
#include <math.h>
#include <estiva/ary.h>
#include <estiva/mesh.h>


static void arrow(FILE *fp,double x0,double y0,double x, double y)
{
  double xl,yl,xr,yr,h;
  h = 0.2;

  xl = (-sqrt(3.0)/2.0)*h; yl = ( 1.0/2.0)*h;
  xr = (-sqrt(3.0)/2.0)*h; yr = (-1.0/2.0)*h;

  fprintf(fp,"%f %f\n",x0,y0);
  fprintf(fp,"%f %f\n",x0+x,y0+y);
  fprintf(fp,"\n");
  fprintf(fp,"%f %f\n",x0+x,y0+y);
  fprintf(fp,"%f %f\n",x0+x+x*xl-y*yl,y0+y+y*xl+x*yl);
  fprintf(fp,"\n");
  fprintf(fp,"%f %f\n",x0+x,y0+y);
  fprintf(fp,"%f %f\n",x0+x+x*xr-y*yr,y0+y+y*xr+x*yr);
  fprintf(fp,"\n");

}




void estiva_pltuv(FILE *fp, xyc *Mid, double *u, double *v)
{
  double t, h, x0,y0,x,y;
  long i,m;
  t = 0.5;
  h = 0.5;
  m = dim1(Mid);
  for(i=1;i<=m;i++){
    x0 = Mid[i].x;
    y0 = Mid[i].y;

    x  = t*u[i];
    y  = t*v[i];

    arrow(fp,x0,y0,x,y);
    fprintf(fp,"%f %f\n",x0,y0);
    fprintf(fp,"%f %f\n",x0+x,y0+y);
    fprintf(fp,"\n");
  }
  fflush(fp);
}

