#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <estiva/op.h>
#include <estiva/ary.h>
#include "msh.h"
#include "spm.h"
#include "date.h"

int incircle(double x, double y, xyc *Z, nde *N,int e2);


#define M( i,j) (*spm_double(M, i,j))
#define K( i,j) (*spm_double(K, i,j))
#define Hx(i,j) (*spm_double(Hx,i,j))
#define Hy(i,j) (*spm_double(Hy,i,j))
#define A( i,j) (*spm_double(A, i,j))

static xyc     *Z, *Mid;
static nde     *N;
static double *S;
static long m, n;

static double length(int a, int b)
{
  double ax, ay, bx, by;
  ax = Z[a].x; ay = Z[a].y;
  bx = Z[b].x; by = Z[b].y;
  return (bx-ax)*(bx-ax)+(by-ay)*(by-ay);
}

static double inner(int a, int b, int c)
{
  double ax, ay, bx, by, cx, cy;
  ax = Z[a].x; ay = Z[a].y;
  bx = Z[b].x; by = Z[b].y;
  cx = Z[c].x; cy = Z[c].y;
  return (cx-ax)*(bx-cx)+(cy-ay)*(by-cy);
}

static spm* M__(void)
{
  static spm *M;
  long  e, a, b, c, A, B, C;

  ary1(M,m+1);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;
    
    M(A,A) += S[e]/3.0;
    M(B,B) += S[e]/3.0;
    M(C,C) += S[e]/3.0;
  }
  return M;
}

static spm* K__(void)
{
  static spm *K;
  long e, a, b, c, A, B, C;
  double s;  

  ary1(K,m+1);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;
    s = S[e];

    K(A,A) += length(b,c)/s; K(A,B) +=inner(a,b,c)/s; K(A,C) +=inner(a,c,b)/s;
    K(B,A) +=inner(b,a,c)/s; K(B,B) += length(c,a)/s; K(B,C) +=inner(b,c,a)/s;
    K(C,A) +=inner(c,a,b)/s; K(C,B) +=inner(c,b,a)/s; K(C,C) += length(a,b)/s;

  }
  return K;
}

static spm* Hx__(void)
{
  static spm *Hx;
  long e, a, b, c, A, B, C;
  
  ary1(Hx,m+1);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;

    Hx(A,e) += -(Z[c].y-Z[b].y);
    Hx(B,e) += -(Z[a].y-Z[c].y);
    Hx(C,e) += -(Z[b].y-Z[a].y);
  }
  return Hx;
}

static spm* Hy__(void)
{
  static spm *Hy;
  long e, a, b, c, A, B, C;
  
  ary1(Hy,m+1);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;

    Hy(A,e) += (Z[c].x-Z[b].x);
    Hy(B,e) += (Z[a].x-Z[c].x);
    Hy(C,e) += (Z[b].x-Z[a].x);

  }
  return Hy;
}

static spm* A__(spm *M, double t, spm *K, spm *Hx, spm *Hy)
{
  static spm *A;
  long   i, j, k, NUM;
  NUM = m*2+n;
  ary1(A, NUM+1);

  for(i=1;i<=m;i++) fornonzero(K[i],j){
    A(  i,   j) = M(i,j) + t*K(i,j);
    A(m+i, m+j) = M(i,j) + t*K(i,j);
  }
  for(i=1;i<=m;i++) fornonzero(Hx[i],j){
    A(i,2*m+j) = -t*Hx(i,j);
    A(2*m+j,i) = -t*Hx(i,j);
  }
  for(i=1;i<=m;i++) fornonzero(Hy[i],j){
    A(m+i,2*m+j) = -t*Hy(i,j);
    A(2*m+j,m+i) = -t*Hy(i,j);
  }
  return A;
}

static double* Fx_(void)
{
  static double *Fx;
  ary1(Fx,m+1);
  return Fx;
}

static double* Fy_(void)
{
  static double *Fy;
  long i;
  ary1(Fy,m+1);
/*
  for(i=1;i<=m;i++) Fy[i] = -9.8; 
*/
  return Fy;
}


static double* Ux_(void)
{
  static double *Ux;
  ary1(Ux,m+1);
  return Ux;
}

static double* Uy_(void)
{
  static double *Uy;
  ary1(Uy,m+1);
  return Uy;
}


static double* b_(spm *M,double t,double *Fx,double *Fy,double *Ux, double *Uy)
{
  static double *b, Mjj;
  long   i, j, k, NUM;
  NUM = m*2+n;
  ary1(b,NUM+1);
  

  for(i=1;i<=m;i++){
    b[i]   = t*M(i,i)*Fx[i] + M(i,i)*Ux[i];
    b[m+i] = t*M(i,i)*Fy[i] + M(i,i)*Uy[i];
  }

  for(i=1;i<=n;i++) b[i+2*m] = 0.0;
  return b;
}

static void boundary_condition(spm *A, double *b)
{
  long NUM;
  int i, j;
  NUM = 2*m+n;

  forgamma(Mid,i,"j1"){
    fornonzero(A[i],j) A(i,j) = 0.0;
    A(i,i) = 1.0;
    b[i] = 0.0;
  }

  forgamma(Mid,i,"j2"){
    fornonzero(A[i],j) A(i,j) = 0.0;
    A(i,i) = 1.0;
    b[i]   = 0.1;

    fornonzero(A[i+m],j) A(i+m,j) = 0.0;
    A(i+m,i+m) = 1.0;
    b[i+m]     = 0.0;
  }

  forgamma(Mid,i,"j3"){
    fornonzero(A[i],j) A(i,j) = 0.0;
    A(i,i) = 1.0;
    b[i]   = 0.0;
  }
  forgamma(Mid,i,"j4"){
    fornonzero(A[i+m],j) A(i+m,j) = 0.0;
    A(i+m,i+m) = 1.0;
    b[i+m]     = 0.0;
  }

  forgamma(Mid,i,"ki"){
    fornonzero(A[i],j) A(i,j) = 0.0;
    A(i,i) = 1.0;
    b[i] = 0.0;

    fornonzero(A[i+m],j) A(i+m,j) = 0.0;
    A(i+m,i+m) = 1.0;
    b[i+m]     = 0.0;
  }

  fornonzero(A[NUM],j) A(NUM,j) = 0.0;
  A(NUM,NUM) = 1.0;
  b[NUM] = 0.0;
/*  ritz(A,b); */
}

static void job(double *p)
{
  FILE *pp;
  long e, a, b, c, i, j;
  pp = stdout;

  
  for(e=1;e<=n;e++){
    a = N[e].a; b = N[e].b; c = N[e].c;
    
    fprintf(pp,"%f %f %f\n",Z[a].x,Z[a].y,p[e]);
    fprintf(pp,"%f %f %f\n",Z[b].x,Z[b].y,p[e]);
    fprintf(pp,"%f %f %f\n",Z[c].x,Z[c].y,p[e]);
    fprintf(pp,"%f %f %f\n",Z[a].x,Z[a].y,p[e]);
    fprintf(pp,"\n\n");
  }
}

void arrow(FILE *fp,double x0,double y0,double x, double y)
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



static void uvplot(FILE *fp,double *u, double *v)
{
  double t, h, x0,y0,x,y,xl,yl,xr,yr;
  long i;
  t = 0.5;
  h = 0.5;
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


static void bplot(FILE *fp, xyc *Z, nde *N,double *P)
{
  long i,a,b,c;
  for(i=1;i<=dim1(N);i++){
    a = N[i].a; b = N[i].b; c = N[i].c;
    if(strcmp(Z[a].label,"inner")&&strcmp(Z[b].label,"inner")){
      fprintf(fp,"%f %f %f\n",Z[a].x,Z[a].y,P[a]);
      fprintf(fp,"%f %f %f\n",Z[b].x,Z[b].y,P[b]);
      fprintf(fp,"\n\n");
    }
    if(strcmp(Z[b].label,"inner")&&strcmp(Z[c].label,"inner")){
      fprintf(fp,"%f %f %f\n",Z[b].x,Z[b].y,P[b]);
      fprintf(fp,"%f %f %f\n",Z[c].x,Z[c].y,P[c]);
      fprintf(fp,"\n\n");
    }
    if(strcmp(Z[c].label,"inner")&&strcmp(Z[a].label,"inner")){
      fprintf(fp,"%f %f %f\n",Z[c].x,Z[c].y,P[c]);
      fprintf(fp,"%f %f %f\n",Z[a].x,Z[a].y,P[a]);
      fprintf(fp,"\n\n");
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

static void p0plot(FILE *fp, xyc *Z, nde *N, double *P)
{
  long e, a, b, c;
  for(e=1;e<=dim1(N);e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    fprintf(fp,"%f %f %f\n",Z[a].x,Z[a].y,P[e]);
    fprintf(fp,"%f %f %f\n",Z[b].x,Z[b].y,P[e]);
    fprintf(fp,"%f %f %f\n",Z[c].x,Z[c].y,P[e]);
    fprintf(fp,"%f %f %f\n",Z[a].x,Z[a].y,P[e]);
    fprintf(fp,"\n\n");
  }
}

static void p1plot(FILE *fp, xyc *Z, nde *N, double *P)
{
  long e, a, b, c;
  for(e=1;e<=dim1(N);e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    fprintf(fp,"%f %f %f\n",Z[a].x,Z[a].y,P[a]);
    fprintf(fp,"%f %f %f\n",Z[b].x,Z[b].y,P[b]);
    fprintf(fp,"%f %f %f\n",Z[c].x,Z[c].y,P[c]);
    fprintf(fp,"%f %f %f\n",Z[a].x,Z[a].y,P[a]);
    fprintf(fp,"\n\n");
  }
}
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
    exit(1);
  }
  for(i=1;i<=dim1(p1);i++) p1[i]/= (double)pi[i];
  return p1;
}

static void frplot(FILE *fp, xyc *Z, nde *N, double scale)
{
  long i, a, b, c;
  double xmax, ymax, xmin, ymin, xmed, ymed, t;
  xmax = Z[1].x, ymax = Z[1].y;
  xmin = Z[1].x, ymin = Z[1].y;
  for(i=1;i<=dim1(N);i++){
    a = N[i].a, b = N[i].b, c = N[i].c;
    
    if(xmax < Z[a].x) xmax = Z[a].x;
    if(xmax < Z[b].x) xmax = Z[b].x;
    if(xmax < Z[c].x) xmax = Z[c].x;
    
    if(xmin > Z[a].x) xmin = Z[a].x;
    if(xmin > Z[b].x) xmin = Z[b].x;
    if(xmin > Z[c].x) xmin = Z[c].x;
    
    if(ymax < Z[a].y) ymax = Z[a].y;
    if(ymax < Z[b].y) ymax = Z[b].y;
    if(ymax < Z[c].y) ymax = Z[c].y;
    
    if(ymin > Z[a].y) ymin = Z[a].y;
    if(ymin > Z[b].y) ymin = Z[b].y;
    if(ymin > Z[c].y) ymin = Z[c].y;
  }
  xmed = ( xmax + xmin )/2.0;
  ymed = ( ymax + ymin )/2.0;
  t = scale;

  fprintf(fp,"%f %f %f\n",xmed-(t*1.3),ymed-t,0.0);
  fprintf(fp,"%f %f %f\n",xmed+(t*1.3),ymed-t,0.0);
  fprintf(fp,"%f %f %f\n",xmed+(t*1.3),ymed+t,0.0);
  fprintf(fp,"%f %f %f\n",xmed-(t*1.3),ymed+t,0.0);
  fprintf(fp,"%f %f %f\n",xmed-(t*1.3),ymed+t,0.0);
}

static void p0top1plot(FILE *fp, xyc *Z, nde *N, double *p)
{
  static double *P;
  P = p0top1(Z,N,p);
  p1plot(fp,Z,N,P);
}

static void pplot(FILE *fp, xyc *Z, nde *N, double *p)
{
  static double *P;
  static long  *Pi;
  double pmax, pmin, med, k, h;
  long i, n, a, b, c;


  P = p0top1(Z,N,p);
  
  pmax = p[1];
  pmin = p[1];
  for(i=2;i<=dim1(N);i++){
    if(pmax < p[i]) pmax = p[i];
    if(pmin > p[i]) pmin = p[i];
  }
  h = 1.0;
  if(defop("-conta"))
    h = atof(getop("-conta"));
  
  med = (pmax + pmin)/2.0;
  fprintf(stderr,"\n(pmax,pmin) = (%f, %f)\n",pmax,pmin);
  conta(fp,P,Z,N,med);
  for(k = med+h; k < pmax; k += h) conta(fp,P,Z,N,k);
  for(k = med-h; k > pmin; k -= h) conta(fp,P,Z,N,k);


  bplot(fp,Z,N,P);
  
  

  fflush(fp);
}

static void pltbound(FILE *fp, xyc *Z, nde *N)
{
  double xo,yo,x,y;
  long e,a,b,c;


  for(e=1;e<=dim1(N);e++){
    a = N[e].a, b = N[e].b, c = N[e].c;

    if(Z[a].label!=NULL)if(Z[b].label!=NULL)if(strcmp(Z[a].label,"inner")&&strcmp(Z[b].label,"inner")){
      xo = Z[a].x; yo = Z[a].y;
      x  = Z[b].x; y  = Z[b].y;
      fprintf(fp,"%f %f\n%f %f\n\n",xo,yo,x,y);
    }
    if(Z[b].label!=NULL)if(Z[c].label!=NULL)if(strcmp(Z[b].label,"inner")&&strcmp(Z[c].label,"inner")){
      xo = Z[b].x; yo = Z[b].y;
      x  = Z[c].x; y  = Z[c].y;
      fprintf(fp,"%f %f\n%f %f\n\n",xo,yo,x,y);
    }
    if(Z[c].label!=NULL)if(Z[a].label!=NULL)if(strcmp(Z[c].label,"inner")&&strcmp(Z[a].label,"inner")){
      xo = Z[c].x; yo = Z[c].y;
      x  = Z[a].x; y  = Z[a].y;
      fprintf(fp,"%f %f\n%f %f\n\n",xo,yo,x,y);
    }
  }
   

  fflush(fp);
}

static void pltdot(FILE *fp, xyc *Z, nde *N)
{
  double xo,yo,x,y, h;
  long e,a,b,c;
  
  h = 0.05;
  for(a=1;a<=dim1(Z);a++)
    fprintf(fp,"%f %f\n%f %f\n\n",
	    Z[a].x,Z[a].y,Z[a].x,Z[a].y);


  fflush(fp);
}

xyc *G_(xyc *Z, nde *N);
int main(int argc, char** argv)
{
  static spm  *M, *K, *Hx, *Hy, *A;
  static double *Fx, *Fy, *Ux, *Uy, *b;
  static xyc * G;
  
  double t=0.1;
  int             i, j, k, e, NUM;
  FILE *fp, *tfp, *pp, *pfp, *ppp;

  static char name[100];

  initop(argc,argv);

  fp = (FILE *)argf(argc,argv);

  fp2msh(fp, &Z, &N);
  fclose(fp);
  pltdot(fopen("pltdot","w"),Z,N);
  pltbound(fopen("pltbound","w"),Z,N);
  pltmsh(fopen("pltmsh","w"),Z,N);
  if(defop("-bound")) exit(1);

  if(defop("-t")) t=atof(getop("-t"));
  fprintf(stderr,"t= %f\n",t);
  Mid = Ver2Mid(Z,N);
  S   = S_(Z,N);
  G   = G_(Z,N);

  m = dim1(Mid);
  n = dim1(S);

  M  = M__();
  K  = K__();
  Hx = Hx__();
  Hy = Hy__();

  Fx = Fx_();
  Fy = Fy_();
  Ux = Ux_();
  Uy = Uy_();

  NUM = 2*m+n;
  
  pp = popen("gnuplot","w");
  ppp= popen("gnuplot","w");
  for(k=1;;k++){
    
    A  = A__(M,t,K,Hx,Hy);
    b = b_(M,t,Fx,Fy,Ux,Uy);
    boundary_condition(A,b);

    ritz(A,b);
    solver(A,b);


    tfp = fopen("flow","w");
    uvplot(tfp,b,&b[m]);
    fclose(tfp);
    fprintf(pp,"plot \"flow\" with lines\n");
    fflush(pp);

    pfp = fopen("pressure","w");
    if(defop("-conta")){
      pplot(pfp,Z,N,&b[2*m]);
      fclose(pfp);
      fprintf(ppp,"plot \"pressure\" with lines\n");
    }
    else{
      p0top1plot(pfp,Z,N,&b[2*m]);
      fclose(pfp);
      fprintf(ppp,"splot \"pressure\" with lines\n");
    }
    fflush(ppp);


    fprintf(stderr,"k = %d",k);
    for(i=1;i<=m;i++){ Ux[i] = b[i]; Uy[i] = b[i+m];}
  }
}


