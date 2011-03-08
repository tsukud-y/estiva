#include "estiva/TaylorHood.h"

#include <math.h>
#include <stdlib.h>

static double B1, B2, C1, C2;
double estiva_TaylorHood_Delta;

void estiva_setBCD(double b1, double b2, double c1, double c2, double s) {
  B1 = b1; B2 = b2; C1 = c1; C2 = c2, Delta = s;
}

double estiva_TaylorHood_a(long i){
  switch(i) {
  case 1: return  0.0;
  case 2: return  0.0;
  case 3: return  1.0;
  case 4: return  0.0;
  case 5: return  0.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}

double estiva_TaylorHood_b(long i){
  switch(i) {
  case 1: return -1.0;
  case 2: return  0.0;
  case 3: return -3.0;
  case 4: return  0.0;
  case 5: return  4.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}

double estiva_TaylorHood_c(long i){
  switch(i) {
  case 1: return  0.0;
  case 2: return -1.0;
  case 3: return -3.0;
  case 4: return  4.0;
  case 5: return  0.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}

double estiva_TaylorHood_d(long i){
  switch(i) {
  case 1: return  0.0;
  case 2: return  0.0;
  case 3: return  4.0;
  case 4: return -4.0;
  case 5: return -4.0;
  case 6: return  4.0;
  default: abort();
  }
  abort();
  return NAN;
}

double estiva_TaylorHood_e(long i){
  switch(i) {
  case 1: return  2.0;
  case 2: return  0.0;
  case 3: return  2.0;
  case 4: return  0.0;
  case 5: return -4.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}

double estiva_TaylorHood_f(long i){
  switch(i) {
  case 1: return  0.0;
  case 2: return  2.0;
  case 3: return  2.0;
  case 4: return -4.0;
  case 5: return  0.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}

double estiva_TaylorHood_alphaB(long j)
{
  switch(j){
  case 1: 
  case 2: 
  case 3: 
  case 4: 
  case 5: 
  case 6: return b(j)*B1 + c(j)*B2;
  default: abort();
  }
  abort();
  return NAN;
}


double estiva_TaylorHood_betaB(long j)
{
  switch(j){
  case 1: 
  case 2: 
  case 3: 
  case 4: 
  case 5: 
  case 6: return 2.0*e(j)*B1 + d(j)*B2;
  default: abort();
  }
  abort();
  return NAN;
}


double estiva_TaylorHood_gammaB(long j)
{
  switch(j){
  case 1: 
  case 2: 
  case 3: 
  case 4: 
  case 5: 
  case 6: return d(j)*B1 + 2.0*f(j)*B2;
  default: abort();
  }
  abort();
  return NAN;
}

double estiva_TaylorHood_alphaC(long j)
{
  switch(j){
  case 1: 
  case 2: 
  case 3: 
  case 4: 
  case 5: 
  case 6: return b(j)*C1 + c(j)*C2;
  default: abort();
  }
  abort();
  return NAN;
}


double estiva_TaylorHood_betaC(long j)
{
  switch(j){
  case 1: 
  case 2: 
  case 3: 
  case 4: 
  case 5: 
  case 6: return 2.0*e(j)*C1 + d(j)*C2;
  default: abort();
  }
  abort();
  return NAN;
}


double estiva_TaylorHood_gammaC(long j)
{
  switch(j){
  case 1: 
  case 2: 
  case 3: 
  case 4: 
  case 5: 
  case 6: return d(j)*C1 + 2.0*f(j)*C2;
  default: abort();
  }
  abort();
  return NAN;
}


double estiva_TaylorHood_ad(long j){
  switch(j) {
  case 1: return  0.0;
  case 2: return  0.0;
  case 3: return  1.0;
  default: abort();
  }
  abort();
  return NAN;
}

double estiva_TaylorHood_bd(long j){
  switch(j) {
  case 1: return  1.0;
  case 2: return  0.0;
  case 3: return -1.0;
  default: abort();
  }
  abort();
  return NAN;
}

double estiva_TaylorHood_cd(long j){
  switch(j) {
  case 1: return  0.0;
  case 2: return  1.0;
  case 3: return -1.0;
  default: abort();
  }
  abort();
  return NAN;
}

