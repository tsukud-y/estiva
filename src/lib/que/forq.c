#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva/que.h"


#define  forall(m,i,n) for(i=m;i<=n;i++) 
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


static que **estiva_p(void **e)
{
  return  f(e);
}


#define p(e)  (*estiva_p(e))

void estiva_forqinit(que *q, void **e)
{
  if (f(e) == NULL) R(e,malloc(sizeof(void *)));
  p(e) = q;
}

int estiva_forq2(void **e)
{
  if ( p(e)->elem == NULL ) { 
    free(f(e)); 
    R(e,NULL);
    return 0; 
  }
  *e = *(void**)p(e)->elem;
  p(e) = p(e)->next;
  return 1;
}
