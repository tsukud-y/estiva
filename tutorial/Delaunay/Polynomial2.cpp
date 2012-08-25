#include "Mesh.h"
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

void Mesh::Polynomial2(vector<Xyc>&Z, vector<Nde>&N)
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

