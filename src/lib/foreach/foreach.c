#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <estiva/foreach.h>

#define  forall(m,i,n) for(i=m;i<=n;i++) 


static void estiva_cp(void *A, void *B, int size1, int size2)
     /* cp(a,b)     estiva_cp(&(a),&(b),sizeof(a),sizeof(b)) */
{ static int i; char *p1, *p2;
  p1=(char *)A;
  p2=(char *)B;

  if(size1 == size2) forall(0,i,size1-1) p2[i]=p1[i]; 
  else{ fprintf(stderr,"estiva_cp: size mismatch\n"); abort();}
}



#define  R(x,y) estiva_R(x,y) 
#define  f(x)   estiva_f(x) 

#define set(i,xi,yi) (x[i]=(xi),y[i]=(yi))

static void *x[20],*y[20];
static int tail=0;

static void *X(int n)
{ return x[n];}

static void *Y(int n)
{ return y[n];}

static int n(void *xi)
{ int i=0; forall(0,i,tail)if(X(i)==xi) break; return i;}

static void estiva_R(void *xi, void *yi)
     /* R(x,y) estiva_R(x,y) */
{ 
  if(yi==NULL)
    { int i;i=n(xi);set(i,NULL,NULL); if(i==tail) tail--;}
  else{
    R(xi,NULL);
    { int i;i=n(NULL);set(i,xi,yi); if(i==tail) tail++;}
  }
}

static void *estiva_f(void *xi)
     /* f(x)   estiva_f(x) */
{ return Y(n(xi));}


#define n(f0) (*(int *)f(f0))

int estiva_foreach(int fsize,void *f0,...)
     /* foreach(f) while(estiva_foreach(sizeof(f),&f, */
{ int i; void *fi, *fn; va_list ap;
  
  if(fsize==0||f0==NULL) return 0;    
  if(f(f0)==NULL){ R(f0,malloc(sizeof(int))); n(f0) = 0;}
  if(f(f0)==NULL) abort();
  
  va_start(ap,f0); fi = f0; 
  forall(1,i,n(f0)) fi = va_arg(ap, void *);
  fn = va_arg(ap,void *);
  va_end(ap); 
  
  estiva_cp(f0,fi,fsize,fsize);
  if(fn!=NULL){ estiva_cp(fn,f0,fsize,fsize); n(f0)++; return 1;}
  free(f(f0)); R(f0,NULL); return 0;
}
void *estiva_end()
     /* end estiva_end())) */
{ return NULL;}
