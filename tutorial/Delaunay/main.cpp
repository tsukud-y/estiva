#include "Delaunay.h"
#undef push
#undef pop
#include <stack>
#include <cmath>
#include <unistd.h>
#include <algorithm>


long FindPair(vector<Nde>&N, long e, long a)
{
  long i, b, c, t, ea, eb, ec;


  ea = N[e].a; eb= N[e].b; ec=N[e].c;
  while( a != ea ) {
    t =ec; ec = eb; eb = ea; ea = t;
  }
  b = eb; c = ec;
  
  for ( i = 1; i < N.size(); i++) {
    if ( e == i ) continue;
    if ( b == N[i].a || b == N[i].b || b == N[i].c ){
      if ( c == N[i].a || c == N[i].b || c == N[i].c )
	return i;
    }
  }
  return 0;
}

void GenRelation(vector<Xyc> Z, vector<Nde>&N)
{
  long i;
  
  for ( i = 1 ; i < N.size(); i++) {
    long a = N[i].a, b = N[i].b, c = N[i].c;
    // printf("Pair is %d %d\n",FindPair(N,i,a),a);
    //printf("Pair is %d %d\n",FindPair(N,i,b),b);
    //printf("Pair is %d %d\n",FindPair(N,i,c),c);

    N[i].A = FindPair(N,i,a);
    N[i].B = FindPair(N,i,b);
    N[i].C = FindPair(N,i,c);
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

  Putmesh(Z,N);
  Xmesh(pp,Z,N);
  sleep(300);
}
