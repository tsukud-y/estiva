#include <stdio.h>
#include "msh.h"

void mshplot(xyc *Z, nde *N)
{
  long i;
  printf("<xyc>\n");
  for(i=1;i<=dim1(Z);i++) printf("%2d    %.14f %.14f %s\n",i,Z[i].x,Z[i].y,Z[i].label);
  printf("<nde>\n");
  for(i=1;i<=dim1(N);i++) printf("%2d    %2d %2d %2d    %2d %2d %2d\n",
				 i,N[i].a,N[i].b,N[i].c,
				 N[i].A,N[i].B,N[i].C);

}
