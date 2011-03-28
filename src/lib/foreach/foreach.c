#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "estiva/foreach.h"
#include "estiva/ary.h"
#include "estiva/std.h"

#define n static_bind(int,x)

int estiva_foreach(int xsize, void *x, ...)
{ 
  va_list ap;
  void *xi, *xn;
  int i; 

  va_start(ap,x); 
  xi = x;
  forall(1,i,n) 
    xi = va_arg(ap, void *);
  xn = va_arg(ap,void *);
  va_end(ap); 
  
  memcpy(xi,x,xsize);
  
  if ( xn != NULL ) { 
    memcpy(x,xn,xsize);
    n++; 
    return 1;
  }
  static_free(x);
  return 0;
}

void *estiva_foreachend()
{ 
  return NULL;
}
