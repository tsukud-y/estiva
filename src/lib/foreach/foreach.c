#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "estiva/foreach.h"
#include "estiva/ary.h"
#include "estiva/std.h"

#define y(x) (*(int *)f(x))

int estiva_foreach(int xsize, void *x, ...)
{ 
  va_list ap;
  void *xi, *xn;
  int i; 
  
  if ( xsize == 0    ) return 0;
  if ( x     == NULL ) return 0;    
  if ( &y(x) == NULL ) { Rnew(x, int);  y(x) = 0; }
  if ( &y(x) == NULL ) abort();

  va_start(ap,x); 
  xi = x;
  forall(1,i,y(x)) 
    xi = va_arg(ap, void *);
  xn = va_arg(ap,void *);
  va_end(ap); 
  
  memcpy(xi,x,xsize);

  if ( xn != NULL ) { 
    memcpy(x,xn,xsize);
    y(x)++; 
    return 1;
  }
  Rdestroy(x);
  return 0;
}

void *estiva_foreachend()
{ 
  return NULL;
}
