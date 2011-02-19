#include <math.h>
#include <stdio.h>
#include <string.h>
#include <estiva/op.h>
#include <estiva/ary.h>
#include <estiva/esolver.h>

/* invert.c : 1999 3 8 by Oguni */
// modified by T. KUDO 
// modified by TSUKUDA Yoshio

static const int Itr=100;
static const double tol=1e-6;

#define Matrix double **
#define Vector double *

static void lu(Matrix A,Matrix L,Matrix U,double tol,int n);
static void subs(Matrix L,Matrix U,Vector b,Vector x,double tol,int n);
static double invert(Matrix A,Vector x,double tol,int n,int *itr);

double estiva_minesolver(double **A, double *x)
{
  int n=dim1(x), itr[2]={Itr, 0};
  return invert(A,x,tol,n,itr);
}


/* function invert */
static double invert(Matrix A,Vector x,double tol,int n,int *itr)
{
  static Matrix L;
  static Matrix U;
  static Vector y;
  static Vector s;
  double c, d, xmin=0.0,lambda=10000000.0;
  int m, i;

  ary2(L,n+1,n+1);
  ary2(U,n+1,n+1);
  ary1(y,n+1);
  ary1(s,n+1);

  lu(A,L,U,tol,n);           // AをLU分解する。結果はLとUに入れる
				  for (m=1; m<=itr[0]; m++){
    for (i=1; i<=n; i++)
       s[i]=x[i]; 
    subs(L,U,s,y,tol,n); //　LUy = s の解yを求める処理。
    c=0; d=0;
    for (i=1; i<=n; i++){
      c+=y[i]*x[i];
      d+=y[i]*y[i];
    }
    lambda =c/d;
    itr[1]=m;
    if (fabs(xmin-lambda)<tol){
      return lambda;  /* 収束した */
    }
    xmin=lambda; 
    for (i=1; i<=n; i++)
        x[i]=y[i]/sqrt(d);  // xにyを正規化して代入
				 }
  itr[1]=m;
  return lambda;  //Itr回では収束しなかった
		      }

/* function lu  (注）Uの対角要素を１とする方法をとっている。*/
static void lu(Matrix A,Matrix L,Matrix U,double tol,int n)
{
  double akk, aik;
  int i, j, k;
  for (k=1; k<=n-1; k++){
    if (fabs(A[k][k])<tol) 
      printf("The pivot is too small\n");
    akk=1.0/A[k][k];
    for (i=k+1; i<=n; i++){
      aik=-A[i][k]*akk;
      for (j=k+1; j<=n; j++)
        A[i][j]+=aik*A[k][j];
    }
    for (j=k+1; j<=n; j++)
      A[k][j]*=akk; 
  }
  for (i=1; i<=n-1; i++){
    for (j=1; j<=i; j++)
      L[i][j]=A[i][j];
    for (j=i+1; j<=n; j++)
      U[i][j]=A[i][j];
    U[i][i]=1.0;
  }
  U[n][n]=1.0;
  for (j=1; j<=n; j++)
    L[n][j]=A[n][j]; 
}
/* function subs */
static void subs(Matrix L,Matrix U,Vector b,Vector x,double tol,int n)
{
  static Vector y;
  double s;
  int i, j, k;

  ary1(y,n+1);

  y[1]=b[1]/L[1][1];
  for (i=2; i<=n; i++){
    if (fabs(L[i][i])<tol) 
      printf("The pivot is too small\n");
    s=0.0;
    for (j=1; j<=i-1; j++)
      s+=L[i][j]*y[j];
    y[i]=(b[i]-s)/L[i][i];
  }
  x[n]=y[n]/U[n][n];
  for (k=n-1; k>=1; k--){
    s=0.0;
    for (j=k+1; j<=n; j++)
      s+=U[k][j]*x[j];
    x[k]=(y[k]-s)/U[k][k];
  }
}
