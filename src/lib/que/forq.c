#include <stdio.h>
#include <stdlib.h>
#include "estiva/que.h"
#include "estiva/std.h"

#define y(x) (*(que**)f(x))

void estiva_forq(que *q, void **x)
{
  if (&y(x) == NULL) Rnew(x,void*);
  if (&y(x) == NULL) abort();
  y(x) = q;
}

int estiva_forq_loop(void **x)
{
  if ( y(x)->elem == NULL ) { 
    Rdestroy(x);
    return 0; 
  }
  *x = *(void**)y(x)->elem;
  y(x) = y(x)->next;
  return 1;
}
