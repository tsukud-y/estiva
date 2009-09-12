#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <estiva/op.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/solver.h>
#include <estiva/mesh.h>


#define M( i,j) mx(M, i,j)
#define K( i,j) mx(K, i,j)
#define Hx(i,j) mx(Hx,i,j)
#define Hy(i,j) mx(Hy,i,j)
#define A( i,j) mx(A, i,j)

static xyc     *Z, *Mid;
static nde     *N;
static double *S;
static long m, n;

extern double *S_(xyc *Z, nde *N);
extern void *Ver2Mid(xyc *Z, nde *N);

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

static MX* M__(void)
{
  static MX *M;
  long  e, a, b, c, A, B, C;

  initmx(M,m+1,20);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;
    
    M(A,A) += S[e]/3.0;
    M(B,B) += S[e]/3.0;
    M(C,C) += S[e]/3.0;
  }
  return M;
}

static MX* K__(void)
{
  static MX *K;
  long e, a, b, c, A, B, C;
  double s;  

  initmx(K,m+1,20);

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

static MX* Hx__(void)
{
  static MX *Hx;
  long e, a, b, c, A, B, C;
  
  initmx(Hx,m+1,20);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;

    Hx(A,e) += -(Z[c].y-Z[b].y);
    Hx(B,e) += -(Z[a].y-Z[c].y);
    Hx(C,e) += -(Z[b].y-Z[a].y);
  }
  return Hx;
}

static MX* Hy__(void)
{
  static MX *Hy;
  long e, a, b, c, A, B, C;
  
  initmx(Hy,m+1,20);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;

    Hy(A,e) += (Z[c].x-Z[b].x);
    Hy(B,e) += (Z[a].x-Z[c].x);
    Hy(C,e) += (Z[b].x-Z[a].x);

  }
  return Hy;
}

static MX* A__(MX *M, double t, MX *K, MX *Hx, MX *Hy)
{
  static MX *A;
  long   i, j, NUM;
  NUM = m*2+n;
  initmx(A, NUM+1,50);

  for(i=1;i<=m;i++) for(j=1; j<=m; j++){
    A(  i,   j) = M(i,j) + t*K(i,j);
    A(m+i, m+j) = M(i,j) + t*K(i,j);
  }
  for(i=1;i<=m;i++) for(j=1; j<=n; j++){
    A(i,2*m+j) = -t*Hx(i,j);
    A(2*m+j,i) = -t*Hx(i,j);
  }
  for(i=1;i<=m;i++) for(j=1; j<=n; j++){
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
  ary1(Fy,m+1);
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


static double* b_(MX *M,double t,double *Fx,double *Fy,double *Ux, double *Uy)
{
  static double *b;
  long   i, NUM;
  NUM = m*2+n;
  ary1(b,NUM+1);
  

  for(i=1;i<=m;i++){
    b[i]   = t*M(i,i)*Fx[i] + M(i,i)*Ux[i];
    b[m+i] = t*M(i,i)*Fy[i] + M(i,i)*Uy[i];
  }

  for(i=1;i<=n;i++) b[i+2*m] = 0.0;
  return b;
}

static void boundary_condition(MX *A, double *b)
{
  long NUM;
  int i, j;
  NUM = 2*m+n;

  forgamma(Mid,i,"j1"){
    for(j=1;j<=NUM;j++) A(i,j) = 0.0;
    A(i,i) = 1.0;
    b[i] = 0.0;
  }

  forgamma(Mid,i,"j2"){
    for(j=1;j<=NUM;j++) A(i,j) = 0.0;
    A(i,i) = 1.0;
    b[i]   = 0.1;

    for(j=1;j<=NUM;j++) A(i+m,j) = 0.0;
    A(i+m,i+m) = 1.0;
    b[i+m]     = 0.0;
  }

  forgamma(Mid,i,"j3"){
    for(j=1;j<=NUM;j++) A(i,j) = 0.0;
    A(i,i) = 1.0;
    b[i]   = 0.0;
  }
  forgamma(Mid,i,"j4"){
    for(j=1;j<=NUM;j++) A(i+m,j) = 0.0;
    A(i+m,i+m) = 1.0;
    b[i+m]     = 0.0;
  }

  forgamma(Mid,i,"ki"){
    for(j=1;j<=NUM;j++) A(i,j) = 0.0;
    A(i,i) = 1.0;
    b[i] = 0.0;

    for(j=1;j<=NUM;j++) A(i+m,j) = 0.0;
    A(i+m,i+m) = 1.0;
    b[i+m]     = 0.0;
  }

  
  for(j=1;j<=NUM;j++) A(NUM,j) = 0.0;
  A(NUM,NUM) = 1.0;
  b[NUM] = 0.0;

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
  double t, h, x0,y0,x,y;
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
  double h;
  long a;
  
  h = 0.05;
  for(a=1;a<=dim1(Z);a++)
    fprintf(fp,"%f %f\n%f %f\n\n",
	    Z[a].x,Z[a].y,Z[a].x,Z[a].y);


  fflush(fp);
}

xyc *G_(xyc *Z, nde *N);
int main(int argc, char** argv)
{
  static MX  *M, *K, *Hx, *Hy, *A;
  static double *Fx, *Fy, *Ux, *Uy, *b;
  static xyc * G;
  
  double t=0.1;
  int             i,  k, NUM;
  FILE *fp, *tfp, *pp, *pfp, *ppp;


  initop(argc,argv);

  fp = (FILE *)stdfp();

  estiva_fp2mesh(fp, &Z, &N);
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


    {
      static double *x;
      static MX *AMX;
      int n, i, j;
      n = dim1(b);
      ary1(x,n+1);
      initmx(AMX,n+1,80);
      for ( i=1; i<=n; i++) for (j=1; j<=n; j++) mx(AMX,i,j) = A(i,j);

      solver(AMX,x,b);
      for ( i=1; i<=n; i++ ) b[i] = x[i];
    }



    tfp = fopen("flow","w");
    uvplot(tfp,b,&b[m]);
    fclose(tfp);
    fprintf(pp,"plot \"flow\" with lines\n");
    fflush(pp);

    pfp = fopen("pressure","w");
    {
      int i;
      static double *p;
      n = dim1(S);
      ary1(p,n+1);
      for ( i = 0; i<=n; i++ ) p[i] = b[2*m+i];
      plt(pfp,Z,N,p);
      fclose(pfp);
      fprintf(ppp,"splot \"pressure\" with lines\n");
    }
    fflush(ppp);


    fprintf(stderr,"k = %d",k);
    for(i=1;i<=m;i++){ Ux[i] = b[i]; Uy[i] = b[i+m];}
  }
}
