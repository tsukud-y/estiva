#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "estiva.h"
#include "Rf.h"
#include "foreach.h"

#define n(f0) (*(int *)f(f0))

int estiva_foreach(int fsize,void *f0,...)
     /* foreach(f) while(estiva_foreach(sizeof(f),&f, */
{ int i; void *fi, *fn; va_list ap;
  
  if(fsize==0||f0==NULL) return 0;    
  if(f(f0)==NULL){ R(f0,malloc(sizeof(int))); n(f0) = 0;}
  if(f(f0)==NULL) exit(1);
  
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

