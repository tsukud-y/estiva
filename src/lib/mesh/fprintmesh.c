#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva/std.h"
#include <estiva/foreach.h>
#include <estiva/ary.h>
#include "estiva/mesh.h"
#include <estiva/op.h>


static char *label(xyc *Z,int i)
{ static char *normal = "inner";
  return Z[i].label==NULL? normal:Z[i].label;
}

void estiva_fprintmesh(void *vfp, xyc *Z, nde *N)
{
  FILE *fp;
  fp = vfp;
  int i;
  fprintf(fp,"<xyc>\n");
  forall(1,i,dim1(Z))
    fprintf(fp,"%4d %.9f %.9f %s\n",i,Z[i].x,Z[i].y,label(Z,i));
  fprintf(fp,"<nde>\n");
  forall(1,i,dim1(N))
    fprintf(fp,"%4d  %4d %4d %4d  %4d %4d %4d\n",
            i, N[i].a,N[i].b,N[i].c,N[i].A,N[i].B,N[i].C);
}

