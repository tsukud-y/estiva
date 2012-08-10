#include "Delaunay.h"
#undef push
#undef pop
#include <stack>
#include <cmath>
#include <unistd.h>

static int isR(Xyc a,Xyc b)
{
  double dot = (a.x*b.x + a.y*b.y);

  if ( fabs(dot) < 0.00001) return 1;
  return 0;
}


long SearchBy2P(long a, long b, vector<Nde>&N)
{
  long i;

  for ( i = 1; i<N.size(); i++) 
    if ( N[i].a == a || N[i].b == a || N[i].c == a )
      if ( N[i].a == b || N[i].b == b || N[i].c == b )
	break;
  return i;
}


void SwapT(long e, long i, vector<Nde>&N)
{
  long a, b, c, x, y, z, t;

  while ( N[e].a == N[i].a || N[e].a == N[i].b || N[e].a == N[i].c ) {
    t =N[e].c; N[e].c = N[e].b; N[e].b = N[e].a; N[e].a = t;
  }

  while ( N[i].a == N[e].a || N[i].a == N[e].b || N[i].a == N[e].c ) {
    t =N[i].c; N[i].c = N[i].b; N[i].b = N[i].a; N[i].a = t;
  }

  a=N[e].a; b=N[e].b; c=N[e].c;
  x=N[i].a; y=N[i].b; z=N[i].c;

  N[e].a = a; N[e].b= x; N[e].c= c;
  N[i].a = x; N[i].b= a; N[i].c= z;
}

void VanishBT(vector<Xyc>&Z,vector<Nde>&N)
{
  for (long  i = 1; i< N.size(); i++ ){
    
    if ( Z[N[i].a].label!=""&&Z[N[i].b].label!=""&&Z[N[i].c].label!= "") {
      Xyc abv, bcv, cav;
      long e;
      abv.x = Z[N[i].b].x-Z[N[i].a].x;  abv.y = Z[N[i].b].y-Z[N[i].a].y; 
      bcv.x = Z[N[i].c].x-Z[N[i].b].x;  bcv.y = Z[N[i].c].y-Z[N[i].b].y; 
      cav.x = Z[N[i].a].x-Z[N[i].c].x;  cav.y = Z[N[i].a].y-Z[N[i].c].y; 

      if (isR(cav,abv)) {
	e = SearchBy2P(N[i].b,N[i].c,N);
      }//A
      if (isR(abv,bcv)) {
	e = SearchBy2P(N[i].c,N[i].a,N);
      }//B
      if (isR(bcv,cav)) {
	e = SearchBy2P(N[i].a,N[i].b,N);
      }//C
      if (incircleDelta(N[e].a,i,Z,N)&&
	  incircleDelta(N[e].b,i,Z,N)&&
	  incircleDelta(N[e].c,i,Z,N)) SwapT(e,i,N);
    }
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

  //Putmesh(Z,N);
  Xmesh(pp,Z,N);
  sleep(300);
}
