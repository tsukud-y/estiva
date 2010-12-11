#include "fem.h"

#include <stdio.h>
#include <string.h>
#include "estiva/que.h"
#include "estiva/mesh.h"
#include "estiva/ary.h"


long *estiva_forgammap2_p, estiva_forgammap2_flag;


que *estiva_forgammap2_init(xyc *Z, nde *N, char *label)
{
  static que *q;
  long i, e, m, *p;

  initq(q);
  forq(q,p) pop(q,p);

  for (i=1; i<=dim1(Z); i++) {
    if (Z[i].label && !strcmp(Z[i].label,label)){
      push(q,i);
    }
  }

  for (e=1; e<=dim1(N); e++) {
    m = N[e].A;
    if ( (Z[N[e].b].label && !strcmp(Z[N[e].b].label, label)) ||
         (Z[N[e].c].label && !strcmp(Z[N[e].c].label, label))  ) {
      if (Z[N[e].b].label && Z[N[e].c].label) push(q,m);
    }
    m = N[e].B;
    if ( (Z[N[e].c].label && !strcmp(Z[N[e].c].label, label)) ||
         (Z[N[e].a].label && !strcmp(Z[N[e].a].label, label))  ) {
      if (Z[N[e].c].label && Z[N[e].a].label) push(q,m);
    }
    m = N[e].C;
    if ( (Z[N[e].a].label && !strcmp(Z[N[e].a].label, label)) ||
         (Z[N[e].b].label && !strcmp(Z[N[e].b].label, label))  ) {
      if (Z[N[e].a].label && Z[N[e].b].label) push(q,m);
    }
  }
  return q;
}