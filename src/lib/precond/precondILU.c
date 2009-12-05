#include <stdio.h>
#include <stdlib.h>
#include <estiva/ary.h>
#include <estiva/solver.h>

#define forall(m,i,n) for(i=m;i<=n;i++)
#define less(x,y)     (x<y?(x):(y))
#define absv(x)       ((x)>0.0?(x):-(x))

static long halfbw(double** a)
{ 
  long n, i, j;
  if(dim2(a) != dim1(a)||dim2(a) < 1||dim1(a) < 1||dim0(a) < 1) return 0;
  
  n = dim2(a);
  for(i=1;i<=n;i++) for(j=1;j<=i;j++)
    if(a[j][n-(i-j)] != 0.0||a[n-(i-j)][j] != 0.0) return n-i;
  return n;
}


static long* LU(double  **a)
{  double nrm,aik,*apk,*ai,apkk;

   static long i,j,k,w,n,*p,pk,pklimit,klimit; 
   n=dim2(a);w=halfbw(a);
   
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


static long gauss_c(long *pivot, double **a, double* b)
{ 
  double *api,s,bpi;
  long *p,i,j,n,w,pi,pilimit; 

  if(dim1(a) != dim1(b)) return 0;

  n = dim2(a);
  w = halfbw(a);

  //if(NULL==(p=LU(a))) return 0;
  //for (i=0; i<=n; i++) printf("p[%d]=%d\n",i,p[i]);
  p = pivot;
  
  for(i=1;i<=n;i++){
    pilimit=less(i+w,n); bpi=b[p[i]];
    forall(i+1,j,pilimit) b[j] += a[j][i]*bpi; 
  }
  for(i=n;i>=1;i--){
    s = b[(pi=p[i])]; api=a[pi]; pilimit=less(pi+w,n);
    forall(i+1,j,pilimit) s -= api[j]*b[j];
    b[pi] = s/api[i];
  }

  {
    static int counter = 0;
    printf("counter = %d\n",counter++);
  }
  return 1;
}

#include <estiva/mx.h>
#define A(i,j) mx(A,i,j)

static long gauss(long *pivot, void* A, double* b)
{
  static double **M;
  long i, j, n;

  n = dim1(b);
  ary2(M,n+1,n+1);

  for(i=1; i<=n; i++) for(j=1; j<=n; j++) M[i][j] = A(i,j);
  return gauss_c(pivot,M,b);
}

long estiva_ILUsolver(long *pivot, void *A, double *x, double *b)
{
  long i, n;
  n = dim1(b);

  for (i=0; i<=n; i++) x[i] = b[i];  
  return gauss(pivot,A,x);
}


void estiva_precondILU(long *pivot, MX *A, double *x, double *b)
{
  long i;
  static double *x1, *b1;

  ary1(x1,A->m+1); ary1(b1,A->m+1);

  x1[0] = 0.0;
  for (i=1; i<=A->m; i++) x1[i] = x[i-1];
    
  b1[0] = 0.0;
  for (i=1; i<=A->m; i++) b1[i] = b[i-1];
  
  estiva_ILUsolver(pivot,A,x1,b1);

  for (i=1; i<=A->m; i++) x[i-1] = x1[i];
  for (i=1; i<=A->m; i++) b[i-1] = b1[i];

  printf("hello precondILU solver\n");
}

void ILU(long *pivot, MX *A)
{
  static long *p;
  static double **ILU;
  long i, j, n;
  
  n = A->m;
  ary2(ILU,n+1,n+1);
  
  for(i=1; i<=n; i++) for(j=1; j<=n; j++) ILU[i][j] = A(i,j);

  p= LU(ILU);

  for(i=0; i<=n; i++) pivot[i] = p[i];

  for(i=1; i<=n; i++) for(j=1; j<=n; j++) A(i,j) = ILU[i][j];

  printf("hello ILU decomposition!V\n");
}
