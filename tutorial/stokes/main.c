#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <estiva/op.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/solver.h>
#include <estiva/mesh.h>
#include <estiva/tmpfile.h>


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
extern void *np1(xyc *Z, nde *N);

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

  fp = stdfp();
  fp2mesh(fp, &Z, &N);
  fclose(fp);

  if(defop("-t")) t=atof(getop("-t"));
  fprintf(stderr,"t= %f\n",t);

  Mid = np1(Z,N);
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
    static double *x, *p;

    A  = A__(M,t,K,Hx,Hy);
    b = b_(M,t,Fx,Fy,Ux,Uy);

    boundary_condition(A,b);

    ary1(x,dim1(b)+1);
    ary1(p,dim1(S)+1);

    solver(A,x,b);
    for (i=1; i<=dim1(p); i++) p[i] = x[2*m+i];

    tfp = tmpopen();
    pltuv(tfp, Mid, x,&x[m]);
    fprintf(pp,"plot \"%s\" w l\n",tmpname(tfp));
    fflush(pp);
    tmpclose(tfp);

    pfp = tmpopen();
    plt(pfp,Z,N,p);
    fprintf(ppp,"splot \"%s\" w l\n",tmpname(pfp));
    fflush(ppp);
    tmpclose(pfp);

    for(i=1;i<=m;i++){ Ux[i] = x[i]; Uy[i] = x[i+m];}
    fprintf(stderr,"k = %d\n",k);
  }
}
