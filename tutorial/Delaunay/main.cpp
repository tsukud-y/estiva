#include "Delaunay.h"
#undef push
#undef pop
#include <stack>
#include <cmath>
#include <unistd.h>
#include <algorithm>



long max(vector<Nde>&N)
{
  long i,max=0;
  for ( i = 1; i<(long)N.size(); i++){
    if ( N[i].a > max ) max = N[i].a;
    if ( N[i].b > max ) max = N[i].b;
    if ( N[i].c > max ) max = N[i].c;;
  }
  return max;
}

void Polynomial2(vector<Xyc>&Z, vector<Nde>&N)
{
  long i, m, n;

  for ( i = 1; i<(long)N.size(); i++) {
    N[i].A *= -1;
    N[i].B *= -1;
    N[i].C *= -1;
  }
  n = max(N)+1;

  long u;
  for(i=1;i<(long)N.size();i++){

    m = N[i].A;
    if(m<=0){ 
      u = N[-m].A;
      if ( u == -i ) u = n;
      N[-m].A = u;
      u = N[-m].B;
      if ( u == -i ) u = n;
      N[-m].B = u;
      u = N[-m].C;
      if ( u == -i ) u = n;
      N[-m].C = u;
      m = n++;
    }
    N[i].A = m;

    m = N[i].B;
    if(m<=0){ 
      u = N[-m].A;
      if ( u == -i ) u = n;
      N[-m].A = u;
      u = N[-m].B;
      if ( u == -i ) u = n;
      N[-m].B = u;
      u = N[-m].C;
      if ( u == -i ) u = n;
      N[-m].C = u;
      m = n++;
    }
    N[i].B = m;

    m = N[i].C;
    if(m<=0){ 
      u = N[-m].A;
      if ( u == -i ) u = n;
      N[-m].A = u;
      u = N[-m].B;
      if ( u == -i ) u = n;
      N[-m].B = u;
      u = N[-m].C;
      if ( u == -i ) u = n;
      N[-m].C = u;
      m = n++;
    }
    N[i].C = m;
  }
}



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

void xmeshp2(vector<Xyc>&Z,vector<Nde>&N) {
if ( fork() == 0 ) {
    long i;
    FILE *pp;
    pp = popen("gnuplot -geometry +0+0","w");
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
    sleep(30*60);
    pclose(pp);
    exit(0);
  }
}


int main(int argc, char ** argv)
{
  vector<Xyc> Z; vector<Nde> N;  
  Xyc z;

  Z.push_back(z);

  for ( z.y = 0.0; z.y <= 1.0; z.y+=0.125)
    for ( z.x = 0.0; z.x <= 1.0; z.x+=0.125)
      {
	if ( z.x == 0.0 || z.y == 0.0 || z.x == 1.0 || z.y == 1.0)
	  z.label="G";
	else z.label="";
	Z.push_back(z);
      }

  FILE *pp = popen("gnuplot","w");

  GenSuperNodes(Z,N);
  DelaunayAlgo(Z,N);

  VanishSuperNodes(Z,N);

  VanishBT(Z,N);

  SortTri(Z,N);

  GenRelation(Z,N);
  Normalization(Z,N);
  Polynomial2(Z,N);

  Putmesh(Z,N);
  //xmeshp2(Z,N);

  Xmesh(pp,Z,N);

  sleep(300);
}

