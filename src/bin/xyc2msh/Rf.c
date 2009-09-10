#include <stdio.h>
#include <stdlib.h>
#include "estiva.h"
#include "Rf.h"
#define set(i,xi,yi) (x[i]=(xi),y[i]=(yi))

static void *x[20],*y[20];
static int tail=0;

static void *X(int n)
{ return x[n];}

static void *Y(int n)
{ return y[n];}

static int n(void *xi)
{ int i=0; forall(0,i,tail)if(X(i)==xi) break; return i;}

void estiva_R(void *xi, void *yi)
     /* R(x,y) estiva_R(x,y) */
{ 
  if(yi==NULL)
    { int i;i=n(xi);set(i,NULL,NULL); if(i==tail) tail--;}
  else{
    R(xi,NULL);
    { int i;i=n(NULL);set(i,xi,yi); if(i==tail) tail++;}
  }
}
void *estiva_f(void *xi)
     /* f(x)   estiva_f(x) */
{ return Y(n(xi));}
