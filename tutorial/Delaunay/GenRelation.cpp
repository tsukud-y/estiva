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
  
  for ( i = 1; i < (long)N.size(); i++) {
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
  
  for ( i = 1 ; i < (long)N.size(); i++) {
    long a = N[i].a, b = N[i].b, c = N[i].c;

    N[i].A = FindPair(N,i,a);
    N[i].B = FindPair(N,i,b);
    N[i].C = FindPair(N,i,c);
  }
}

