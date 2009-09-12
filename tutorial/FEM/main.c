#include <stdio.h>
#include <math.h>
#include <string.h>
#include <estiva/op.h>
#include <estiva/ary.h>
#include <estiva/mesh.h>
#include <estiva/esolver.h>



static void
set_A(xyc *Z, nde *N, double **A)
{
  int n, i, j, k;
  double S, D, xi, yi, xj, yj, xk, yk, b1, b2, b3, c1, c2, c3;

  for(n=1;n<=dim1(N);n++){
    i=N[n].a; j=N[n].b; k=N[n].c;

    xi = Z[i].x; yi = Z[i].y;
    xj = Z[j].x; yj = Z[j].y;
    xk = Z[k].x; yk = Z[k].y;

    D = xi*(yj-yk)+xj*(yk-yi)+xk*(yi-yj);
    S = fabs(D)/2.0;

    b1=(yj-yk)/D; b2=(yk-yi)/D; b3=(yi-yj)/D;
    c1=(xk-xj)/D; c2=(xi-xk)/D; c3=(xj-xi)/D;

    A[i][i]+=S*(b1*b1+c1*c1);A[i][j]+=S*(b2*b1+c2*c1);A[i][k]+=S*(b3*b1+c3*c1);
    A[j][i]+=S*(b1*b2+c1*c2);A[j][j]+=S*(b2*b2+c2*c2);A[j][k]+=S*(b3*b2+c3*c2);
    A[k][i]+=S*(b1*b3+c1*c3);A[k][j]+=S*(b2*b3+c2*c3);A[k][k]+=S*(b3*b3+c3*c3);
  }

  for(i=1;i<=dim1(Z);i++){
    if(!strcmp("boundary",Z[i].label)) A[i][i]= 1000000000000000000000000000.0;
  }
}



static void
set_u(xyc *Z, double *u)
{
  int i;
  for(i=1; i<=dim1(Z); i++) u[i] = 1.0;
}



main(int argc, char **argv){
  static xyc *Z;
  static nde *N;
  static double **A, *u;

  initop(argc, argv);
  fp2mesh(stdfp(),&Z, &N);

  ary2(A,dim1(Z)+1, dim1(Z)+1); ary1(u,dim1(Z)+1);

  set_A(Z,N,A); set_u(Z,u);

  esolver(A,u);

  plt(NULL,Z,N,u); sleep(1000);
}
