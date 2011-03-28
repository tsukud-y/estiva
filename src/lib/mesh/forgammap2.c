#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva/mesh.h"
#include "estiva/ary.h"
#include "estiva/std.h"
#include "estiva/que.h"

#define y(x) (*(que**)f(x))

void estiva_forgammap2(long *x, xyc *Z, nde *N, char *label)
{
  long v, e, m;
  if (&y(x) == NULL) Rnew(x,void*);
  if (&y(x) == NULL) abort();
  initq(y(x));

  for (v=1; v<=dim1(Z); v++) {
    if (Z[v].label && !strcmp(Z[v].label,label)){
      push(y(x),v);
    }
  }

  for (e=1; e<=dim1(N); e++) {
    m = N[e].A;
    if ( (Z[N[e].b].label && !strcmp(Z[N[e].b].label, label)) ||
         (Z[N[e].c].label && !strcmp(Z[N[e].c].label, label))  ) {
      if (Z[N[e].b].label && Z[N[e].c].label) push(y(x),m);
    }
    m = N[e].B;
    if ( (Z[N[e].c].label && !strcmp(Z[N[e].c].label, label)) ||
         (Z[N[e].a].label && !strcmp(Z[N[e].a].label, label))  ) {
      if (Z[N[e].c].label && Z[N[e].a].label) push(y(x),m);
    }
    m = N[e].C;
    if ( (Z[N[e].a].label && !strcmp(Z[N[e].a].label, label)) ||
         (Z[N[e].b].label && !strcmp(Z[N[e].b].label, label))  ) {
      if (Z[N[e].a].label && Z[N[e].b].label) push(y(x),m);
    }
  }
}

int estiva_forgammap2_loop(long *x)
{
  if( y(x)->elem ) { 
    pop(y(x),*x);
    return 1;
  }
  Rdestroy(x);
  return 0;
}
