#include "Mesh.h"
#undef push
#undef pop
#include <stack>
#include <cmath>
#include <unistd.h>
#include <algorithm>



static void pltmsh(FILE *fp, vector<Xyc>&Z,vector<Nde>&N)
{
  long e, a, b, c;
  for(e=1;e<(long)N.size();e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    fprintf(fp,"%f %f\n",Z[a].x,Z[a].y);
    fprintf(fp,"%f %f\n",Z[b].x,Z[b].y);
    fprintf(fp,"%f %f\n",Z[c].x,Z[c].y);
    fprintf(fp,"%f %f\n",Z[a].x,Z[a].y);
    fprintf(fp,"\n\n");
  }


}

void Mesh::XP2(FILE *pp, vector<Xyc>&Z,vector<Nde>&N) {
if ( fork() == 0 ) {
    long i;

    {
      fprintf(pp,"set border 0\n");
      fprintf(pp,"set format x \"\"\n");
      fprintf(pp,"set format y \"\"\n");
      fprintf(pp,"set noxtics \n");
      fprintf(pp,"set noytics \n");
    }
    {
      for ( i=1; i<(long)Z.size(); i++)
        fprintf(pp,"set label \"%ld\" at %f,%f\n",i,Z[i].x,Z[i].y);
    }
    {
      for ( i=1; i<(long)N.size(); i++)
        fprintf(pp,"set label \"%ld\" at %f,%f\n", N[i].C,
                (Z[N[i].a].x + Z[N[i].b].x)/2.0,
                (Z[N[i].a].y + Z[N[i].b].y)/2.0);
      for ( i=1; i<(long)N.size(); i++)
        fprintf(pp,"set label \"%ld\" at %f,%f\n", N[i].A,
                (Z[N[i].b].x + Z[N[i].c].x)/2.0,
                (Z[N[i].b].y + Z[N[i].c].y)/2.0);
      for ( i=1; i<(long)N.size(); i++)
        fprintf(pp,"set label \"%ld\" at %f,%f\n", N[i].B,
                (Z[N[i].c].x + Z[N[i].a].x)/2.0,
                (Z[N[i].c].y + Z[N[i].a].y)/2.0);

    }
    {
      for ( i=1; i<(long)N.size(); i++)
        fprintf(pp,"set label \"(%ld)\" at %f,%f\n",i,
                (Z[N[i].a].x + Z[N[i].b].x + Z[N[i].c].x)/3.0,
                (Z[N[i].a].y + Z[N[i].b].y + Z[N[i].c].y)/3.0);
    }
    fprintf(pp,"plot '-' title \"\" with lines\n");
    pltmsh(pp,Z,N);
    fprintf(pp,"e\n");
    fflush(pp);
  }
}

