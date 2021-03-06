#include "estivaplus/Mesh.h"                                                                                

void Mesh::X(FILE *fp, vector<Xyc> &Z, vector<Nde> &N)
{
  if ( defop("-XP2") ) {
    Mesh::XP2(fp,Z,N);
    return;
  }

  long e, a, b, c;

  fprintf(fp,"unset label\n");

  for (e=1; e<(long)N.size(); e++) {
    a = N[e].a, b = N[e].b, c = N[e].c;
    fprintf(fp,"set label \"%ld%s\" at %f , %f\n",
	    a,Z[a].label.c_str(),Z[a].x,Z[a].y);
    fprintf(fp,"set label \"%ld%s\" at %f , %f\n",
	    b,Z[b].label.c_str(),Z[b].x,Z[b].y);
    fprintf(fp,"set label \"%ld%s\" at %f , %f\n",
	    c,Z[c].label.c_str(),Z[c].x,Z[c].y);
    fprintf(fp,"set label \"(%ld)\" at %f , %f\n",
	    e,(Z[a].x+Z[b].x+Z[c].x)/3.0,(Z[a].y+Z[b].y+Z[c].y)/3.0);
  }

  fprintf(fp,"plot '-' title \"\" with lines\n");

  for(e=1;e<(long)N.size();e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    fprintf(fp,"%f %f\n",Z[a].x,Z[a].y);
    fprintf(fp,"%f %f\n",Z[b].x,Z[b].y);
    fprintf(fp,"%f %f\n",Z[c].x,Z[c].y);
    fprintf(fp,"%f %f\n",Z[a].x,Z[a].y);
    fprintf(fp,"\n\n");
  }
  fprintf(fp,"e\n");
  fflush(fp);
}
