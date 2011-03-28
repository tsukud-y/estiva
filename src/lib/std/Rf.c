#include <stdio.h>
#include <stdlib.h>
#include "estiva/ary.h"
#include "estiva/std.h"

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

void estiva_std_R(void *x, void *y)
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

void *estiva_std_f(void *x)
{ 
  ary1(x_array,limit); ary1(y_array,limit);
  return y_array[where(x)];
}

void estiva_std_Rnew(void *x,size_t size)
{
  R(x,calloc(1,size));
}

void estiva_std_Rdestroy(void *x)
{
  free(f(x));
  R(x, NULL);
}

