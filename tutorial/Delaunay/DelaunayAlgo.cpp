#include "Mesh.h"
#undef push
#undef pop
#include <stack>


void Mesh::DelaunayAlgo(vector<Xyc>&Z,vector<Nde>&N)
{
  for (long i = 1 ; i <(long)Z.size()-3; i++) {
    stack<long> st;
    long e0,e1,e2;
    static FILE *pp=NULL;
    if (pp == NULL) pp = popen("gnuplot","w");

    e0 = Mesh::SearchT(Z,N,i);
    Mesh::SplitT(i,e0,N);

    if ( 0 != N[e0].A ) st.push(N[e0].A);
    if ( 0 != N[e0].B ) st.push(N[e0].B);
    if ( 0 != N[e0].C ) st.push(N[e0].C);

    st.push(e0);
  
    while(!st.empty()){
      if ( defop("-XMeshAnime") ) Mesh::X(pp,Z,N);


      for ( e1 =0; e1 == 0; ) { e1=st.top(); st.pop();}

      if ( N[e1].a == i ) e2 = N[e1].A;
      if ( N[e1].b == i ) e2 = N[e1].B;
      if ( N[e1].c == i ) e2 = N[e1].C;
      
      if (e2 != 0 )
	if(Mesh::incircle(i,e2,Z,N))
	  if(!Mesh::degeneracy(e1,e2,Z,N)){ st.push(e1); st.push(e2); }
    }
  }
}

