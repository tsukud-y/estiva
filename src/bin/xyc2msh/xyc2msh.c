#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva/std.h"
#include "stack.h"
#include <estiva/foreach.h>
#include <estiva/ary.h>
#include "Delaunay.h"
#include "FILE.h"

#include <estiva/op.h>
static char *label(xyc *Z,int i)
{ static char *normal = "inner";
  return Z[i].label==NULL? normal:Z[i].label;
}
int main(int argc, char **argv)
{ FILE *fp; xyc *Z; nde *N; int i;
  initop(argc,argv);

  if(argc == 1) fp = stdin;
  else if(NULL==(fp=fopen(argv[1],"r"))){
    fprintf(stderr,"%s - version 0.11\n",argv[0]);
    fprintf(stderr,"Can't open file %s\n",argv[1]);
    fprintf(stderr,"usage: %% %s [filename]\n",argv[0]);
    exit(1);
  }
  
  Z = fp2Z(fp); 
  fclose(fp);

  Delaunay(Z,N);

  fp=stdout;
  fprintf(fp,"<xyc>\n");
  forall(1,i,dim1(Z))
    fprintf(fp,"%4d %.9f %.9f %s\n",i,Z[i].x,Z[i].y,label(Z,i));
  fprintf(fp,"<nde>\n");
  forall(1,i,dim1(N))
    fprintf(fp,"%4d  %4d %4d %4d  %4d %4d %4d\n",
	    i, N[i].a,N[i].b,N[i].c,N[i].A,N[i].B,N[i].C);
  fclose(fp);
  return 0;
}
