#include "Delaunay.h"
#undef push
#undef pop
#include <stack>
#include <cmath>
#include <unistd.h>
#include <algorithm>


long min3(long a, long b, long c)
{
  if ( a < b && a < c ) return a;
  if ( b < a && b < c ) return b;
  return c;
}

void Normalization(vector<Xyc>&Z,vector<Nde>&N)
{
  long min, a, b, c, A, B, C, t, i;
  
  for ( i = 1; i < (long)N.size(); i++) {

    a = N[i].a; b = N[i].b; c= N[i].c;  A = N[i].A; B = N[i].B; C = N[i].C;
    
    min = min3(a,b,c);
    while ( min != a ) {
      t = c; c = b; b = a; a = t;
      t = C; C = B; B = A; A = t;
    }
    
    N[i].a = a; N[i].b = b; N[i].c = c; N[i].A = A; N[i].B = B; N[i].C = C;
  }
}


