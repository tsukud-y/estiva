#include <stdio.h>
#include <stdlib.h>
#include "estiva/que.h"
#include "estiva/std.h"

#define p static_bind(que*,x)

void estiva_forq(que *q, void **x)
{
  p = q;
}

int estiva_forq_loop(void **x)
{
  if ( p->elem == NULL ) { 
    static_free(x);
    return 0; 
  }
  *x = *(void**)p->elem;
  p = p->next;
  return 1;
}
