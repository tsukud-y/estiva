#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "estiva/op.h"
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/solver.h"
#include "estiva/mesh.h"
#include "estiva/tmpfile.h"


#define length(a,b) \
((Z[b].x-Z[a].x)*(Z[b].x-Z[a].x)+(Z[b].y-Z[a].y)*(Z[b].y-Z[a].y))

#define inner(a,b,c) \
((Z[c].x-Z[a].x)*(Z[b].x-Z[c].x)+(Z[c].y-Z[a].y)*(Z[b].y-Z[c].y))

static double *S_(xyc *Z, nde *N)
{
  static double *S;
  long i, e;
  double xi,xj,xk,yi,yj,yk;

  e = dim1(N);
  ary1(S,e+1);

  for(i=1;i<=e;i++){
    xi = Z[N[i].a].x;
    xj = Z[N[i].b].x;
    xk = Z[N[i].c].x;
    yi = Z[N[i].a].y;
    yj = Z[N[i].b].y;
    yk = Z[N[i].c].y;
    S[i] = xi*yj+xj*yk+xk*yi-yi*xj-yj*xk-yk*xi;
    S[i] /= 2.0;
  }
  return S;
}

static MX* M__(xyc *Mid,nde *N, double *S)
{
  static MX *M;
  long  e, a, b, c, A, B, C, m, n;

  m = dim1(Mid);
  n = dim1(N);
  initmx(M,m+1,20);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;
    
    mx(M,A,A) += S[e]/3.0;
    mx(M,B,B) += S[e]/3.0;
    mx(M,C,C) += S[e]/3.0;
  }
  return M;
}

static MX* K__(xyc *Mid, xyc *Z, nde *N, double *S)
{
  static MX *K;
  long e, a, b, c, A, B, C, m, n;
  double s;  

  m = dim1(Mid);
  n = dim1(N);
  initmx(K,m+1,20);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;
    s = S[e];

    mx(K,A,A) += length(b,c)/s; mx(K,A,B) +=inner(a,b,c)/s; mx(K,A,C) +=inner(a,c,b)/s;
    mx(K,B,A) +=inner(b,a,c)/s; mx(K,B,B) += length(c,a)/s; mx(K,B,C) +=inner(b,c,a)/s;
    mx(K,C,A) +=inner(c,a,b)/s; mx(K,C,B) +=inner(c,b,a)/s; mx(K,C,C) += length(a,b)/s;

  }
  return K;
}

static MX* Hx__(xyc *Mid, xyc *Z, nde *N)
{
  static MX *Hx;
  long e, a, b, c, A, B, C, m, n;

  m = dim1(Mid);
  n = dim1(N);
  initmx(Hx,m+1,20);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;

    mx(Hx,A,e) += -(Z[c].y-Z[b].y);
    mx(Hx,B,e) += -(Z[a].y-Z[c].y);
    mx(Hx,C,e) += -(Z[b].y-Z[a].y);
  }
  return Hx;
}

static MX* Hy__(xyc *Mid, xyc *Z, nde *N)
{
  static MX *Hy;
  long e, a, b, c, A, B, C, m, n;

  m = dim1(Mid);
  n = dim1(N);
  initmx(Hy,m+1,20);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;

    mx(Hy,A,e) += (Z[c].x-Z[b].x);
    mx(Hy,B,e) += (Z[a].x-Z[c].x);
    mx(Hy,C,e) += (Z[b].x-Z[a].x);

  }
  return Hy;
}

static MX* A__(xyc *Mid, nde *N, MX *M, double t, MX *K, MX *Hx, MX *Hy)
{
  static MX *A;
  long   i, j, NUM, m, n;
  m = dim1(Mid);
  n = dim1(N);
  NUM = m*2+n;
  initmx(A, NUM+1,8);

  for(i=1;i<=m;i++) for(j=1; j<=m; j++){
    mx(A,  i,   j) = mx(M,i,j) + t*mx(K,i,j);
    mx(A,m+i, m+j) = mx(M,i,j) + t*mx(K,i,j);
  }
  for(i=1;i<=m;i++) for(j=1; j<=n; j++){
    mx(A,i,2*m+j) = -t*mx(Hx,i,j);
    mx(A,2*m+j,i) = -t*mx(Hx,i,j);
  }
  for(i=1;i<=m;i++) for(j=1; j<=n; j++){
    mx(A,m+i,2*m+j) = -t*mx(Hy,i,j);
    mx(A,2*m+j,m+i) = -t*mx(Hy,i,j);
  }
  return A;
}


static double* b_(xyc *Mid, nde *N,MX *M,double t,double *Fx,double *Fy,double *Ux, double *Uy)
{
  static double *b;
  long   i, NUM, m, n;

  m = dim1(Mid);
  n = dim1(N);

  NUM = m*2+n;
  ary1(b,NUM+1);
  

  for(i=1;i<=m;i++){
    b[i]   = t*mx(M,i,i)*Fx[i] + mx(M,i,i)*Ux[i];
    b[m+i] = t*mx(M,i,i)*Fy[i] + mx(M,i,i)*Uy[i];
  }

  for(i=1;i<=n;i++) b[i+2*m] = 0.0;
  return b;
}

static void boundary_condition(nde *N, xyc *Mid, MX *A, double *b)
{
  long NUM;
  int i, j, m, n;
  m = dim1(Mid);
  n = dim1(N);
  NUM = 2*m+n;

  forgamma(Mid,i,"j1"){
    for(j=1;j<=NUM;j++) mx(A,i,j) = 0.0;
    mx(A,i,i) = 1.0;
    b[i] = 0.0;
  }

  forgamma(Mid,i,"j2"){
    for(j=1;j<=NUM;j++) mx(A,i,j) = 0.0;
    mx(A,i,i) = 1.0;
    b[i]   = 0.1;

    for(j=1;j<=NUM;j++) mx(A,i+m,j) = 0.0;
    mx(A,i+m,i+m) = 1.0;
    b[i+m]     = 0.0;
  }

  forgamma(Mid,i,"j3"){
    for(j=1;j<=NUM;j++) mx(A,i,j) = 0.0;
    mx(A,i,i) = 1.0;
    b[i]   = 0.0;
  }
  forgamma(Mid,i,"j4"){
    for(j=1;j<=NUM;j++) mx(A,i+m,j) = 0.0;
    mx(A,i+m,i+m) = 1.0;
    b[i+m]     = 0.0;
  }

  forgamma(Mid,i,"ki"){
    for(j=1;j<=NUM;j++) mx(A,i,j) = 0.0;
    mx(A,i,i) = 1.0;
    b[i] = 0.0;

    for(j=1;j<=NUM;j++) mx(A,i+m,j) = 0.0;
    mx(A,i+m,i+m) = 1.0;
    b[i+m]     = 0.0;
  }

  
  for(j=1;j<=NUM;j++) mx(A,NUM,j) = 0.0;
  mx(A,NUM,NUM) = 1.0;
  b[NUM] = 0.0;

}


int main(int argc, char** argv)
{
  static MX  *M, *K, *Hx, *Hy, *A;
  static double *Fx, *Fy, *Ux, *Uy, *b, *x, *p, *S;
  static xyc *Z, *Mid; static nde *N;
  
  double t=0.1;
  int             i,  k, NUM, m, n;
  FILE *fp, *tfp, *pp, *pfp, *ppp;

  initop(argc,argv);

  fp = stdfp();
  fp2mesh(fp, Z, N);
  fclose(fp);

  if(defop("-t")) t=atof(getop("-t"));
  fprintf(stderr,"t= %f\n",t);

  Mid = np1(Z,N);
  S   = S_(Z,N);

  m = dim1(Mid);
  n = dim1(S);
  NUM = 2*m+n;

  M  = M__(Mid,N,S);
  K  = K__(Mid,Z,N,S);
  Hx = Hx__(Mid,Z,N);
  Hy = Hy__(Mid,Z,N);

  ary1(Fx,m+1);  ary1(Fy,m+1);  ary1(Ux,m+1);  ary1(Uy,m+1);
  ary1(x,NUM+1); ary1(p,n+1);
  
  pp = popen("gnuplot","w");
  ppp= popen("gnuplot","w");


  for(k=1;;k++){

    A = A__(Mid, N, M,t,K,Hx,Hy);
    b = b_(Mid,N,M,t,Fx,Fy,Ux,Uy);

    boundary_condition(N,Mid,A,b);

    solver(A,x,b);
    for (i=1; i<=dim1(p); i++) p[i] = x[2*m+i];

    tfp = tmpopen();
    plt(tfp,Mid,Z,N,x);
    fprintf(pp,"plot \"%s\" w l\n",tmpname(tfp));
    fflush(pp);
    tmpclose(tfp);

    pfp = tmpopen();
    plt(pfp,Mid,Z,N,p);
    fprintf(ppp,"splot \"%s\" w l\n",tmpname(pfp));
    fflush(ppp);
    tmpclose(pfp);

    for(i=1;i<=m;i++){ Ux[i] = x[i]; Uy[i] = x[i+m];}
    fprintf(stderr,"k = %d\n",k);
  }
}
