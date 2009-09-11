#include <stdio.h>
#include <math.h>
#include <string.h>
#include <estiva/op.h>
#include <estiva/ary.h>
#include <estiva/mesh.h>
#include <estiva/esolver.h>


void estiva_plt(FILE *fp, xyc *Z, nde *N, double *u)
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
      sleep(1000);
    }
  }
}

static void
set_A(xyc *Z, nde *N, double **A)
{
  int n, i, j, k;
  double S, D, x1, y1, x2, y2, x3, y3;

  double Aii, Aij, Aik;
  double Aji, Ajj, Ajk;
  double Aki, Akj, Akk;

  double b1, b2, b3;
  double c1, c2, c3;

  for(n=1;n<=dim1(N);n++){
    x1 = Z[N[n].a].x; y1 = Z[N[n].a].y;
    x2 = Z[N[n].b].x; y2 = Z[N[n].b].y;
    x3 = Z[N[n].c].x; y3 = Z[N[n].c].y;

    D = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2);
    S = fabs(D)/2.0;

    b1=(y2-y3)/D; b2=(y3-y1)/D; b3=(y1-y2)/D;
    c1=(x3-x2)/D; c2=(x1-x3)/D; c3=(x2-x1)/D;
    
    Aii=S*(b1*b1+c1*c1); Aij=S*(b2*b1+c2*c1); Aik=S*(b3*b1+c3*c1);
    Aji=S*(b1*b2+c1*c2); Ajj=S*(b2*b2+c2*c2); Ajk=S*(b3*b2+c3*c2); 
    Aki=S*(b1*b3+c1*c3); Akj=S*(b2*b3+c2*c3); Akk=S*(b3*b3+c3*c3);
    
    i=N[n].a; j=N[n].b; k=N[n].c;

    A[i][i]+=Aii; A[i][j]+=Aij; A[i][k]+=Aik;
    A[j][i]+=Aji; A[j][j]+=Ajj; A[j][k]+=Ajk;
    A[k][i]+=Aki; A[k][j]+=Akj; A[k][k]+=Akk;
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
  static double **A, *u, lambda;

  initop(argc, argv);
  fp2mesh(stdfp(),&Z, &N);

  ary2(A,dim1(Z)+1, dim1(Z)+1); ary1(u,dim1(Z)+1);

  set_A(Z,N,A); set_u(Z,u);
  
  lambda = esolver(A,u);
  fprintf(stderr,"labmda=%f\n",lambda);

  estiva_plt(NULL,Z,N,u); 
}
