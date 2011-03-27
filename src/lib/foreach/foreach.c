#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "estiva/foreach.h"
#include "estiva/std.h"
#include "estiva/ary.h"

static int top, limit = 1;
static void **x_array, **y_array;

static void set(int i, void *x, void *y)
{
  x_array[i] = x; y_array[i] = y;
}

static int where(void *x)
{ 
  int i; 
  forall(0,i,top) if( x_array[i] == x ) break; 
  return i;
}

static void R(void *x, void *y)
{ 
  int i;

  if ( y == NULL ) { 
    i = where(x);
    set(i,NULL,NULL); 
    if ( i == top ) top--;
  }
  else{
    R(x,NULL);
    i = where(NULL);
    set(i,x,y); 
    if ( i == top ) {
      top++;
      if ( top >= limit ) {
	static void **tx, **ty;
	
	ary1(tx,limit); ary1(ty,limit);
	
	forall(0,i,top) tx[i] = x_array[i];
	forall(0,i,top) ty[i] = y_array[i];

	limit++;
	
	ary1(x_array,limit); ary1(y_array,limit);

	forall(0,i,top) x_array[i] = tx[i];
	forall(0,i,top) y_array[i] = ty[i];
      }
    }
  }
}

static void *f(void *x)
{ 
  ary1(x_array,limit); ary1(y_array,limit);
  return y_array[where(x)];
}

#define y(x) (*(int *)f(x))

static void new(void *x)
{
  R(x,malloc(sizeof(int)));
  y(x) = 0;
}

static void destroy(void *x)
{
  free(&y(x));
  R(x, NULL);
}

int estiva_foreach(int xsize, void *x, ...)
{ 
  va_list ap;
  void *xi, *xn;
  int i; 
  
  if ( xsize == 0    ) return 0;
  if ( x     == NULL ) return 0;    
  if ( &y(x) == NULL ) new(x);
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
  destroy(x);
  return 0;
}

void *estiva_foreachend()
{ 
  return NULL;
}
