#include <stdio.h>
#include "estiva/mesh.h"
#include "estiva/op.h"


int main(int argc, char **argv)
{ 
  FILE *fp; xyc *Z; nde *N; 
  
  initop(argc,argv);
  
  fp = stdfp();
  fp2xyc(fp,Z); 
  fclose(fp);
  
  delaunay(Z,N);
  
  fp = ofp();
  fprintmesh(fp,Z,N);
  fclose(fp);
  return 0;
}
