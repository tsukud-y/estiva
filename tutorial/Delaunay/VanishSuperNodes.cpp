#include "Delaunay.h"

void VanishSuperNodes(vector<Xyc>&Z,vector<Nde>&N)
{
  long e, n=Z.size()-4;


  for ( e = 1; e<N.size(); e++) 
    if (N[e].a == 0 || N[e].b == 0 || N[e].c == 0 ||
	N[e].a >  n || N[e].b >  n || N[e].c >  n ) {
      N[e].a = 0; N[e].b = 0; N[e].c = 0; N[e].A = 0; N[e].B = 0; N[e].C = 0;
    }
  Z.pop_back();
  Z.pop_back();
  Z.pop_back();
  Z[0].x = 0.0; Z[0].y = 0.0; Z[0].label="";

  
  long i, j;
  j = 1;
  for ( i = 1; i+j< N.size();) {
    if (N[i].a==0&& N[i].b ==0 && N[i].c == 0 &&
	N[i].A==0&& N[i].B ==0 && N[i].C == 0 ){
      static Nde zero;
      N[i] = N[i+j]; N[i+j] = zero; j++;
    }
    else { i++; j=0; }
  } 

  for ( j=N.size(); i+1<j; i++)
    N.pop_back();

  for ( i = 0; i<N.size(); i++)
    N[i].A = N[i].B = N[i].C = 0;
}

