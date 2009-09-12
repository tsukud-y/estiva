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
  double S, D, x1, y1, x2, y2, x3, y3, b1, b2, b3, c1, c2, c3;

  for(n=1;n<=dim1(N);n++){
    x1 = Z[N[n].a].x; y1 = Z[N[n].a].y;
    x2 = Z[N[n].b].x; y2 = Z[N[n].b].y;
    x3 = Z[N[n].c].x; y3 = Z[N[n].c].y;

    D = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2);
    S = fabs(D)/2.0;

    b1=(y2-y3)/D; b2=(y3-y1)/D; b3=(y1-y2)/D;
    c1=(x3-x2)/D; c2=(x1-x3)/D; c3=(x2-x1)/D;
    
    i=N[n].a; j=N[n].b; k=N[n].c;

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
