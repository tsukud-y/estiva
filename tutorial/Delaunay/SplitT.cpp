#include "Mesh.h"

void Mesh::SplitT(long i, long e0, vector<Nde>&N)
{
  long e1, e2;
  e1 = N.size();
  e2 = N.size()+1;

  Nde nde;

  N.push_back(nde);
  N.push_back(nde);


  long a,b,c,A,B,C;


  a=N[e0].a; b=N[e0].b; c=N[e0].c; A=N[e0].A; B=N[e0].B; C=N[e0].C;   
  

#define N_set(n,i,j,k,I,J,K) \
  N[n].a=i;N[n].b=j;N[n].c=k;N[n].A=I;N[n].B=J;N[n].C=K;

  N_set( 0,0,0,0, 0, 0, 0);
  N_set(e0,i,b,c,A ,e1,e2);  
  N_set(e1,a,i,c,e0,B, e2);
  N_set(e2,a,b,i,e0,e1, C);

  if(N[B].A==e0) N[B].A=e1;
  if(N[B].B==e0) N[B].B=e1;
  if(N[B].C==e0) N[B].C=e1;

  if(N[C].A==e0) N[C].A=e2;
  if(N[C].B==e0) N[C].B=e2;
  if(N[C].C==e0) N[C].C=e2;
}

