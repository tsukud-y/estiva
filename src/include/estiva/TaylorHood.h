#ifndef _ESTIVA_TAYLORHOOD_H_
#define _ESTIVA_TAYLORHOOD_H_

extern double estiva_TaylorHood_Delta;

double estiva_TaylorHood_a(long i);
double estiva_TaylorHood_b(long i);
double estiva_TaylorHood_c(long i);
double estiva_TaylorHood_d(long i);
double estiva_TaylorHood_e(long i);
double estiva_TaylorHood_f(long i);
double estiva_TaylorHood_ad(long j);
double estiva_TaylorHood_bd(long j);
double estiva_TaylorHood_cd(long j);
double estiva_TaylorHood_alphaB(long j);
double estiva_TaylorHood_betaB(long j);
double estiva_TaylorHood_gammaB(long j);
double estiva_TaylorHood_alphaC(long j);
double estiva_TaylorHood_betaC(long j);
double estiva_TaylorHood_gammaC(long j);

#define a(i)                  estiva_TaylorHood_a(i) 
#define b(i)                  estiva_TaylorHood_b(i) 
#define c(i)                  estiva_TaylorHood_c(i) 
#define d(i)                  estiva_TaylorHood_d(i) 
#define e(i)                  estiva_TaylorHood_e(i) 
#define f(i)                  estiva_TaylorHood_f(i) 
#define ad(j)                 estiva_TaylorHood_ad(j)
#define bd(j)                 estiva_TaylorHood_bd(j)
#define cd(j)                 estiva_TaylorHood_cd(j)
#define alphaB(j)             estiva_TaylorHood_alphaB(j)
#define  betaB(j)             estiva_TaylorHood_betaB(j)
#define gammaB(j)             estiva_TaylorHood_gammaB(j)
#define alphaC(j)             estiva_TaylorHood_alphaC(j)
#define  betaC(j)             estiva_TaylorHood_betaC(j)
#define gammaC(j)             estiva_TaylorHood_gammaC(j)
#define Delta                 estiva_TaylorHood_Delta

#endif
