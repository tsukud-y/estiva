#include <stdio.h>
#include <stdlib.h>
#include <estiva/ary.h>
#include <estiva/mx.h>

#define forall(m,i,n) for(i=m;i<=n;i++)
#define less(x,y)     (x<y?(x):(y))
#define absv(x)       ((x)>0.0?(x):-(x))

#define A(i,j) mx(A,i,j)
static void *A;

static long halfbw(void)
{ 
  long n, i, j;
  
  n = dim1(A);
  for(i=1;i<=n;i++) for(j=1;j<=i;j++)
    if(A(j,n-(i-j)) != 0.0||A(n-(i-j),j) != 0.0) return n-i;
}


static long* LU(double  **a)
{  double nrm,aik,*apk,*ai,apkk;

   static long i,j,k,w,n,*p,pk,pklimit,klimit; 
   n=dim1(A);w=halfbw();
   
   ary1(p,n+1);

   for(k=1;k<=n;k++){
     pk=k; nrm=absv(a[k][k]);
     
     forall(k+1,i,less(k+w,n))if(nrm<absv(a[i][k])){
       pk=i; nrm=absv(a[i][k]);
     }
     if(nrm == 0.0) return NULL;
     p[k]=pk;
     klimit=less(k+w,n); apkk= -a[pk][k];
     forall(k+1,i,klimit){
       a[i][k] /= apkk; 
       pklimit=less(pk+w,n); aik=a[i][k]; apk=a[pk]; ai=a[i];
       forall(k+1,j,pklimit) ai[j] += aik*apk[j];
     }
   }
   return p;
}


static long gauss_c(double **a, double* b)
{ 
  double *api,s,**x,bpi;
  long *p,i,j,k,n,w,pi,pilimit; 

  if(dim1(A) != dim1(b)) return 0;

  n = dim1(A);
  w = halfbw();

  if(NULL==(p=LU(a))) return 0;
  
  for(i=1;i<=n;i++){
    pilimit=less(i+w,n); bpi=b[p[i]];
    forall(i+1,j,pilimit) b[j] += a[j][i]*bpi; 
  }
  for(i=n;i>=1;i--){
    s = b[(pi=p[i])]; pilimit=less(pi+w,n);
    forall(i+1,j,pilimit) s -= a[pi][j]*b[j];
    b[pi] = s/a[pi][i];
  }
  
  return 1;
}



static long blu(void* Ap, double* b)
{
  static double **M, *W;
  long i, j, n, w;
  A = Ap;

  n = dim1(b);
  ary2(M,n+1,1);
  
  w = halfbw();
  ary1(W,(n+2)*(2*w+1));
  
  for(i=0;i<=n;i++) M[i] = &W[i*2*w];

  for(i=1; i<=n; i++) for(j=i-w; j<=i+w; j++) if(1<=j&&j<=n) M[i][j] = A(i,j);
  return gauss_c(M,b);
}


long estiva_blusolver(void *A, double *x, double *b)
{
  long i, n;
  n = dim1(b);

  for (i=0; i<=n; i++) x[i] = b[i];
  return blu(A,x);
}
