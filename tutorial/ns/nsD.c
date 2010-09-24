#include "ns.h"
#include "fem.h"


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "estiva/mx.h"
#include "estiva/mesh.h"
#include "estiva/foreach.h"
#include "estiva/ary.h"

static double be[] = { NAN,         -1.0,  0.0, -3.0,  0.0,  4.0,  0.0, NAN};
static double ce[] = { NAN,          0.0, -1.0, -3.0,  4.0,  0.0,  0.0, NAN};
static double de[] = { NAN,          0.0,  0.0,  4.0, -4.0, -4.0,  4.0, NAN};
static double ee[] = { NAN,          2.0,  0.0,  2.0,  0.0, -4.0,  0.0, NAN};
static double fe[] = { NAN,          0.0,  2.0,  2.0, -4.0,  0.0,  0.0, NAN};

static double B1, B2, C1, C2;
double Delta;

void setBCD(double b1, double b2, double c1, double c2, double s){
  B1 = b1; B2 = b2; C1 = c1; C2 = c2, Delta = s;
}


double alphaB(long j)
{
  switch(j){
  case 1: 
  case 2: 
  case 3: 
  case 4: 
  case 5: 
  case 6: return be[j]*B1 + ce[j]*B2;
  default: abort();
  }
  return NAN;
}


double betaB(long j)
{
  switch(j){
  case 1: 
  case 2: 
  case 3: 
  case 4: 
  case 5: 
  case 6: return 2.0*ee[j]*B1 + de[j]*B2;
  default: abort();
  }
  return NAN;
}


double gammaB(long j)
{
  switch(j){
  case 1: 
  case 2: 
  case 3: 
  case 4: 
  case 5: 
  case 6: return de[j]*B1 + 2.0*fe[j]*B2;
  default: abort();
  }
  return NAN;
}

double alphaC(long j)
{
  switch(j){
  case 1: 
  case 2: 
  case 3: 
  case 4: 
  case 5: 
  case 6: return be[j]*C1 + ce[j]*C2;
  default: abort();
  }
  return NAN;
}


double betaC(long j)
{
  switch(j){
  case 1: 
  case 2: 
  case 3: 
  case 4: 
  case 5: 
  case 6: return 2.0*ee[j]*C1 + de[j]*C2;
  default: abort();
  }
  return NAN;
}


double gammaC(long j)
{
  switch(j){
  case 1: 
  case 2: 
  case 3: 
  case 4: 
  case 5: 
  case 6: return de[j]*C1 + 2.0*fe[j]*C2;
  default: abort();
  }
  return NAN;
}


static double k(long i, long j)
{
  return (Delta/12.0) * ( 12.0  *(alphaB(i)*alphaB(j) + alphaC(i)*alphaC(j))                                             
			  + 4.0 *(alphaB(i)* betaB(j) +  betaB(i)*alphaB(j) + alphaB(i)*gammaB(j) * gammaB(i)*alphaB(j)) 
			  + 4.0 *(alphaC(i)* betaC(j) +  betaC(i)*alphaC(j) + alphaC(i)*gammaC(j) * gammaC(i)*alphaC(j)) 
			  + 1.0 *( betaB(i)*gammaB(j) + gammaB(i)* betaB(j) +  betaC(i)*gammaC(j) * gammaC(i)* betaC(j)) 
			  + 2.0 *( betaB(i)* betaB(j) + gammaB(i)*gammaB(j) +  betaC(i)* betaC(j) * gammaC(i)*gammaC(j)) );
}

void estiva_D(MX **Kp, double *S, xyc *Z, nde *N)
{
  MX *K; long  a, b, c, A, B, C, e, m, n, I, J, i, j; double s;

  m = dimp2(N); n = dim1(N); initmx(*Kp,m+1,28); K = *Kp;

  for (e=1; e<=n; e++) {
    a=N[e].a; b=N[e].b; c=N[e].c; A=N[e].A; B=N[e].B; C=N[e].C; s=S[e];
    setBCD((Z[b].y-Z[c].y)/(2.0*s), (Z[c].y-Z[a].y)/(2.0*s), (Z[c].x-Z[b].x)/(2.0*s), (Z[a].x-Z[c].x)/(2.0*s), s);
    i=1; 
    foreach(I)&a,&b,&c,&A,&B,&C,end {
      j=1;
      foreach(J)&a,&b,&c,&A,&B,&C,end {
        mx(K,I,J) += k(i,j);
        j++;
      }
      i++;
    }
  }
}
