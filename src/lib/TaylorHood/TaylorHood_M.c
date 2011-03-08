#include "estiva/TaylorHood.h"

static double mij(long i, long j)
{
  return  (Delta/180.0)*
    (
     180.0 * (  a(i)*a(j)                                                                                       )
     +60.0  * (  a(i)*b(j) + b(i)*a(j) + a(i)*c(j) + c(i)*a(j)                                                  )
     +15.0  * (  a(i)*d(j) + d(i)*a(j) + b(i)*c(j) + c(i)*b(j)                                                  )
     +30.0  * (  b(i)*b(j) + c(i)*c(j) + a(i)*e(j) + e(i)*a(j) + a(i)*f(j) + f(i)*a(j)                          )
     +18.0  * (  b(i)*e(j) + e(i)*b(j) + c(i)*f(j) + f(i)*c(j)                                                  )
     + 6.0  * (  b(i)*d(j) + d(i)*b(j) + c(i)*e(j) + e(i)*c(j) + b(i)*f(j) + f(i)*b(j) + c(i)*d(j) + d(i)*c(j)  )
     + 2.0  * (  d(i)*d(j) + e(i)*f(j) + f(i)*e(j)                                                              )
     + 3.0  * (  d(i)*e(j) + e(i)*d(j) + d(i)*f(j) + f(i)*d(j)                                                  )
     +12.0  * (  e(i)*e(j) + f(i)*f(j)                                                                          )
     );
}

void estiva_TaylorHood_M(MX **M, long w) {
  genP2P2mx(*M,mij,w);
}
