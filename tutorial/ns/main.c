#include "ns.h"

void   estiva_AX(MX **AXp, double *U, double *S, xyc *Z, nde *N);
void   estiva_AY(MX **AYp, double *V, double *S, xyc *Z, nde *N);
#define nsAX(AX,U,S,Z,N)      estiva_AX(&AX,U,S,Z,N)
#define nsAY(AY,V,S,Z,N)      estiva_AY(&AY,V,S,Z,N)


int main(int argc, char **argv)
{
  static xyc *Z; static nde *N; 
  static MX *A, *K, *M, *Hx, *Hy, *Ax, *Ay; static double *x, *b;
  long   k, kn = 1, m, n, NUM;
  double t = 0.001;

  initop(argc,argv);
  rectmesh(Z,N);

  m   = dimp2(N); 
  n   = dim1(Z); 
  NUM = m*2 + n;

  ary1(x,NUM+1); 
  ary1(b,NUM+1);

  setmesh(Z,N);

  TaylorHood_M(M,16);
  TaylorHood_K(K,12);
  TaylorHood_Hx(Hx,5);
  TaylorHood_Hy(Hy,5);

  if ( defop("-kn") ) kn = atoi(getop("-kn"));
  for ( k = 1; k <= kn; k++ ) {
    TaylorHood_Ax(Ax,x,24);
    TaylorHood_Ay(Ay,x+m,24);
    printf("A\n");
    nsA(A,x,b,K,M,Hx,Hy,Ax,Ay,t,50);
    nsRhs(b,M,x);
    boundary_condition(A,b);
    printf("solver\n");
    solver(A,x,b);
    pltp2(x,Z,N);
    printf("kn - k = %ld\n",kn-k);
   }
  return 0;
}
