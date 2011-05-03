#include "ns.h"


int main(int argc, char **argv)
{
  static xyc *Z; static nde *N; 
  static MX *A, *K, *M, *Hx, *Hy, *Ax, *Ay; static double *x, *b;
  long   k, kn = 100, m, n, NUM;

  initop(argc,argv);
  cylindermesh(&Z,&N);
  fprintmesh(stdout,Z,N);

  m   = dimp2(N); 
  n   = dim1(Z); 
  NUM = m*2 + n;

  ary1(x,NUM+1); 
  ary1(b,NUM+1);

  setmesh(Z,N);


  TaylorHood_M(M,36);

  TaylorHood_K(K,32);

  TaylorHood_Hx(Hx,15);
  TaylorHood_Hy(Hy,15);

  if ( defop("-kn") ) kn = atoi(getop("-kn"));

  for ( k = 1; k <= kn; k++ ) {
    TaylorHood_Ax(Ax,x,34);
    TaylorHood_Ay(Ay,x+m,34);
    nsA(A,x,b,K,M,Hx,Hy,Ax,Ay,50);
    nsRhs(b,M,x);
    boundary_condition(A,b);
    printf("hello1\n");
    solver(A,x,b);
    printf("hello\n");
    pltp2(x,Z,N);
    //estiva_thinplt(x,Z,N);
   }
  return 0;
}
