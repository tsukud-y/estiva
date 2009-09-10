#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include "estiva.h"
#include "confary.h"

#define CENTURY 2001
 
static int history=1, sup=1999;
static void *x[CENTURY], *y[CENTURY];

static void R(void *xi, void *yi)
{ static int i;

  if(yi != NULL){
    R(xi,(void *)NULL);
    forall(1,i,history)if(x[i] == NULL) break;
    x[i] = xi;
    y[i] = yi;
    if(i == history) history++;
    if(sup < history) (void)fprintf(stderr,"The world is over.\n"); 
  }
  else{
    forall(1,i,history)if(x[i] == xi) break;
    x[i] = NULL;
    y[i] = NULL;
  }
}


struct rcs { int d2, d1, sz;};

static struct rcs *f(void *xi) 
{ static int i;
  forall(1,i,history)if(x[i] == xi) break;
  return (struct rcs *)(y[i]);
}

void *confary_ary2(void *q, int m, int n, int size)
{ static int i; void ***p; void **A;
  cp(q,p);cp(*p,A);

  if(dim2(A)==m-1&&dim1(A)==n-1) return A;
  if(A != NULL){ 
    free(f(A)); R(A,NULL); free(A); *p=(void **)NULL;
  }
  if(m < 1 || n < 1 || size < 1) return NULL;

  *p=(void **)calloc(m,sizeof(void *));
  cp(*p,A);
  if(A==NULL) return NULL;
  R(A,malloc((size_t)sizeof(struct rcs)));

  f(A)->d2=m-1, f(A)->d1=n-1, f(A)->sz=size;
  A[0]=calloc(m*n,size);
  return A;
}

void *confary_ary1(void *A, int n, int size)
{ void **p;
  cp(A,p);

  if(dim1(*p) == n-1) return *p;
  if(*p != NULL){ free(f(*p)); R(*p,NULL); free((char *)*p); *p=NULL;}
  if(n < 1 || size < 1) return NULL;
  
  if(NULL == (*p = (void *)calloc((size_t)n,(size_t)size)))
    { fprintf(stderr,"Can't malloc() \n\n");return NULL;}
  
  R(*p,malloc((size_t)sizeof(struct rcs)));
  f(*p)->d2 = -1, f(*p)->d1 = n-1, f(*p)->sz = size;
  reset(*p);
  return *p;
}

int confary_dim2(void *p)
{ if((void *)f(p) != NULL) return f(p)->d2; return -1;}

int confary_dim1(void *p)
{ if((void *)f(p) != NULL) return f(p)->d1; return -1;}

int confary_siz(void *p)
{ if((void *)f(p) != NULL) return f(p)->sz; return 0;}

void confary_reset(void *A)
{ void **p;
  cp(A,p);

  if(siz(p) < 1) return;
  if(dim2(p) < 0) memset((char *)p,'\0',(dim1(p)+1)*siz(p));
  else memset((char *)p[0],'\0',(dim2(p)+1)*(dim1(p)+1)*siz(p));
}

