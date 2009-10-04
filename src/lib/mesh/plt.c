#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <estiva/op.h>
#include <estiva/ary.h>
#include <estiva/mesh.h>

static double *p0top1(xyc *Z, nde *N,double *p0)
{
  static double *p1;
  static long   *pi;
  long i, a, b, c;

  ary1(p1, dim1(Z)+1);
  ary1(pi, dim1(Z)+1);

  for(i=1;i<=dim1(p1);i++){ p1[i] = 0.0; pi[i] = 0;}

  for(i=1;i<=dim1(N);i++){
    a = N[i].a; b = N[i].b; c = N[i].c;

    p1[a] += p0[i];  pi[a]++;
    p1[b] += p0[i];  pi[b]++;
    p1[c] += p0[i];  pi[c]++;
  }

  for(i=1;i<=dim1(p1);i++) if(pi[i]==0){
    fprintf(stderr,"mesh data error\n");
    abort();
  }
  for(i=1;i<=dim1(p1);i++) p1[i]/= (double)pi[i];
  return p1;
}


static void estiva_p1pltinternal(FILE *fp, xyc *Z, nde *N, double *u)
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
  fflush(fp);
}

static void estiva_p1plt(FILE *fp, xyc *Z, nde *N, double *u)
{
  if ( fp == NULL ) {
    fp = fopen("/tmp/plt.tmp","w");
    estiva_p1pltinternal(fp, Z, N, u);
    fclose(fp);
    fp = popen("gnuplot","w");
    fprintf(fp,"splot '/tmp/plt.tmp' w l\n");
    fflush(fp);
  } else  estiva_p1pltinternal(fp, Z, N, u);
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
      int i, n;

      for (i=2; i<=dim1(Z); i++) {
	if(umax < u[i]) umax = u[i];
	if(umin > u[i]) umin = u[i];
      }
      h = (umax-umin)/20.0;

      if( 0.0 != atof(getop("-conta")) ) 
	h = (umax-umin)/(1.0+atof(getop("-conta")));

      n = 1 + atoi(getop("-conta"));
      for (x = umin, i=0; i<=n; x+=h, i++) conta(fp,u,Z,N,x);

      fclose(fp);
      fp = popen("gnuplot","w");
      fprintf(fp,"splot '/tmp/plt.tmp' w l\n");
      fflush(fp);
    }
  } else {
      double x, h, umax = u[1], umin = u[1];
      int i, n;

      for (i=2; i<=dim1(Z); i++) {
	if(umax < u[i]) umax = u[i];
	if(umin > u[i]) umin = u[i];
      }
      h = (umax-umin)/20.0;

      if( 0.0 != atof(getop("-conta")) ) 
	h = (umax-umin)/(1.0+atof(getop("-conta")));
      

      n = 1 + atoi(getop("-conta"));
      for (x = umin, i=0; i<=n; x+=h, i++) conta(fp,u,Z,N,x);

      fflush(fp);
  }
}


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



static void estiva_pltuv(FILE *fp, xyc *Mid, double *u, double *v)
{
  double t, x0,y0,x,y;
  long i,m;
  t = 0.5;

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

void estiva_plt(FILE *fp, xyc *Mid, xyc *Z, nde *N, double *u)
{
  if (defop("-conta")) {
  
    if ( dim1(N) == dim1(u) )
      estiva_contaplt(fp,Z,N,p0top1(Z,N,u));
    else
      estiva_contaplt(fp,Z,N,u);
    return ;
  }
  
  if ( dim1(N) == dim1(u) )
    estiva_p1plt(fp,Z,N,p0top1(Z,N,u));
  else if ( dim1(Z) == dim1(u) )
    estiva_p1plt(fp,Z,N,u);
  else 
    estiva_pltuv(fp,Mid,u,&u[dim1(Mid)]);
}
