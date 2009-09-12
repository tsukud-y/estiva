#include <stdio.h>
#include <estiva/ary.h>
#include <estiva/mesh.h>
#include "forgamma.h"


static int count=0;

void numgrid_initgamma(void)
{ count = 0;}

static int neq(char *str1, char *str2)
{
  if(str1 == (char *)NULL && str2 == (char *)NULL) return 0;
  if(str1 == (char *)NULL || str2 == (char *)NULL) return 1;
  return strcmp(str1,str2);
}

int numgrid_gamma(xyc* Z,int *p,char *str)
{ 
  if(dim1(Z) < ++count) return 0;
  
  while(neq(Z[count].label,str))if(count<dim1(Z)) count++;else return 0;
  *p = count;
  return 1;
}
