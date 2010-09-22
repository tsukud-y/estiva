#include "fem.h"


#include "estiva/ary.h"
#include "estiva/std.h"


long estiva_dimp2(nde *N)
{
  long e, E, n=1;
  E = dim1(N);
  for ( e = 1; e <= E; e++ ) {
    n = max(n,N[e].A);
    n = max(n,N[e].B);
    n = max(n,N[e].C);
  }
  return n;
}
