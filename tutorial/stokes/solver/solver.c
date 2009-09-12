#include <stdio.h>
#include <string.h>
#include <estiva/op.h>
#include "ary.h"
#include "solver.h"

void solver(void *A, double *b)
{
  char *method;

  if(defop("-solver")){
    method = getop("-solver");
    if(!strcmp(method,"gauss")) gauss(A,b);
    if(!strcmp(method,"pcg")) pcg(A,b);
    if(!strcmp(method,"blu")) blu(A,b);
    if(!strcmp(method,"pcgs")) pcgs(A,b);
    if(!strcmp(method,"etalu")) etalu(A,b);
    if(!strcmp(method,"pcgslow")) pcgslow(A,b);
  }
  else gauss(A,b);
}
