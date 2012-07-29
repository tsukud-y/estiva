#include "Delaunay.h"
#undef push
#undef pop
#include <stack>
#include <cmath>
#include <unistd.h>


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
}


int main(int argc, char ** argv)
{
  vector<Xyc> Z; vector<Nde> N;  
  Xyc z;

  z.x=0.0, z.y=0.0, z.label="super_node"; Z.push_back(z);

  for ( z.y = 0.0; z.y < 1.0; z.y+=0.1)
    for ( z.x = 0.0; z.x < 1.0; z.x+=0.1)
      {
	z.label="";
	Z.push_back(z);
      }

  GenSuperNodes(Z,N);

  FILE *pp = popen("gnuplot","w");

  for (long i = 1 ; i <Z.size()-3; i++) {
    stack<long> st;
    long e0,e1,e2;

    e0 = SearchT(Z,N,i);
    SplitT(i,e0,N);

    if ( 0 != N[e0].A ) st.push(N[e0].A);
    if ( 0 != N[e0].B ) st.push(N[e0].B);
    if ( 0 != N[e0].C ) st.push(N[e0].C);

    st.push(e0);
  
    while(!st.empty()){

      for ( e1 =0; e1 == 0; ) { e1=st.top(); st.pop();}

      if ( N[e1].a == i ) e2 = N[e1].A;
      if ( N[e1].b == i ) e2 = N[e1].B;
      if ( N[e1].c == i ) e2 = N[e1].C;
      
      if (e2 != 0 )
	if(incircle(i,e2,Z,N))
	  if(!degeneracy(e1,e2,Z,N)){ st.push(e1); st.push(e2); }
    }
  }

  VanishSuperNodes(Z,N);

  Putmesh(Z,N);
  Xmesh(pp,Z,N);
  sleep(30);
}
