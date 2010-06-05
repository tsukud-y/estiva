#include <stdio.h>
#include <unistd.h>
#include "estiva/op.h"
#include "estiva/mesh.h"
#include "estiva/ary.h"
#include <unistd.h>
#include <sys/types.h>


int main(int argc, char **argv) 
{
  FILE *fp;
  xyc  *Z;

  initop(argc,argv);
  fp = stdfp();
  fp2xyc(fp,Z);
  fclose(fp);
  
  xmesh(Z);
  
  return 0;
}
