#include <stdio.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/op.h"
#include "estiva/solver.h"
#include "estiva/std.h"
#include "estiva/vec.h"
#include "estiva/eblas.h"


int estiva_symcheckmx(void *Apointer)
{
  MX *A;
  long i, j;
  A = Apointer;

  fornonzeromx (A,i,j) if( mx(A,i,j) != mx(A,j,i) ) return 0;
  return 1;
}
