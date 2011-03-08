#include "estiva/TaylorHood.h"

static double hyij(long i, long j)
{
  return (Delta/12.0) * (12.0*(alphaC(i)*ad(j))                                                       
			 +4.0*( betaC(i)*ad(j) + alphaC(i)*bd(j) + gammaC(i)*ad(j) + alphaC(i)*cd(j)) 
			 +1.0*( betaC(i)*cd(j) + gammaC(i)*bd(j))                                     
			 +2.0*( betaC(i)*bd(j) + gammaC(i)*cd(j))                                     );
}

void estiva_TaylorHood_Hy(MX **Hy, long w) {
  genP2P1mx(*Hy,hyij, w);
}
