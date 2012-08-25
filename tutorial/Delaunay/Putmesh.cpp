#include "Delaunay.h"

void fputMesh(FILE *fp,vector<Xyc> &Z, vector<Nde> &N)
{
  fprintf(fp,"<xyc>\n");
  for (unsigned long i=1; i<Z.size(); i++)
    fprintf(fp,"%f %f %s\n",Z[i].x,Z[i].y,Z[i].label.c_str());

  fprintf(fp,"<nde>\n");
  for (unsigned long i=1; i<N.size(); i++)
    fprintf(fp,"%ld %ld %ld %ld %ld %ld\n",
           N[i].a,N[i].b,N[i].c,N[i].A,N[i].B,N[i].C);

 
}
