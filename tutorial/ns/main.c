#include <stdio.h>
#include <stdlib.h>

#include "estiva/ary.h"
#include "estiva/op.h"
#include "estiva/solver.h"

#include "fem.h"
#include "ns.h"



int main(int argc, char **argv)
{
  static xyc *Z; static nde *N; 
  static MX *A, *K, *M, *Hx, *Hy, *AX, *AY; static double *x, *b, *S;
  long   k, kn = 1, m, n, NUM;
  double t = 0.001;

  initop(argc,argv);
  rectmesh(Z,N);

  m   = dimp2(N); 
  n   = dim1(Z); 
  NUM = m*2 + n;

  ary1(x,NUM+1); 
  ary1(b,NUM+1);

  femdelta(S,Z,N); 
  setZNS(Z,N,S);
  genP2P2mx(M,mij);
  genP2P2mx(K,kij);
  genP2P1mx(Hx,hxij);
  genP2P1mx(Hy,hyij);

  
  if ( defop("-kn") ) kn = atoi(getop("-kn"));
  for ( k = 1; k <= kn; k++ ) {
    nsAX(AX,x,S,Z,N);
    nsAY(AY,x+m,S,Z,N);
    nsA(A,x,b,Z,N,K,M,Hx,Hy,AX,AY,t);
    nsRhs(b,M,x);
    boundary_condition(Z,N,A,b);
    solver(A,x,b);
    clearmx(AX);
    clearmx(AY);
    printf("kn - k = %ld\n",kn-k);
    pltp2(x,Z,N);
   }
  return 0;
}
