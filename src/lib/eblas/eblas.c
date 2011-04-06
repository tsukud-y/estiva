#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/std.h"
#include "estiva/vec.h"
#include "estiva/op.h"
#include "estiva/eblas.h"
#undef min
#undef max

static MX *A, *AT;

void *estiva_setAmx(void *Apointer)
{
  A = Apointer;
  return A;
}

void *estiva_setATmx(void *ATpointer)
{
  AT = ATpointer;
  return AT;
}

int estiva_matvec(double *alpha, double *x, double *beta, double *y)
{
  matvecmx(A, alpha, x, beta, y); 
  return 0;
}

int estiva_matvectrans(double *alpha, double *x, double *beta, double *y)
{
  matvecmx(AT, alpha, x, beta, y); 
  return 0;
}

int estiva_psolveq(double *x, double *b, char *str, int L)
{ 
  cpvec(b,x);
  if (!strcmp(str,"LEFT") ) {
    return 0;
  }
  psolvevec(A,x);
  return 0;
}

int estiva_psolvetransq(double *x, double *b, char *str, int L)
{ 
  cpvec(b,x);
  if (!strcmp(str,"LEFT") ) {
    return 0;
  }
  psolvevec(AT,x);
  return 0;
}

int estiva_daxpy_(int *n, double *alpha, double *x, int *c1, double *y, int *c2)
{
  setveclength(*n);
  addvec(*alpha,x,y);
  return 0;
}

int estiva_dcopy_(int *n,double  *b, int *c__1, double *r, int *c__2)
{
  setveclength(*n);
  cpvec(b,r);
  return 0;
}

double estiva_ddot_(int *n, double *z, int *c__1, double *y, int *c__2)
{
  setveclength(*n);
  return dotvec(z,y);
} 

double estiva_dnrm2_(int *n, double *x, int *c__1)
{
  setveclength(*n);
  return L2(x);
}


int estiva_dscal_(int *n, double *d__1, double *v, int *c__1)
{
  long i;
  setveclength(*n);
  forall(0,i,*n-1) v[i] *= *d__1;
  return 0;
}

double estiva_getbreak_(void)
{
  return 1.2e-31;
}

int estiva_psolve(double *x, double *b)
{
  setveclength(A->n);
  cpvec(b,x);
  psolvevec(A,x);
  return 0;
}

static double psc98(MX *A, double *x, double *b)
{
  static double *L;
  static long itr = 0;
  long i;
  ary1(L,A->n+1);
  matvecvec(A,1.0,x,0.0,L);

  if ( defop("-redview") ) {
    forall(0,i,A->n) printf("%ld %ld %e\n",itr,i,fabs(b[i]-L[i]));
    itr++;
    printf("\n");
  }
  return Linf(addformula( L, '=', L, '-',1.0,b));
}

static double B = 0.0;

double estiva_setpsc98Linf(double Bvalue)
{
  B = Bvalue;
  return B;
}

int estiva_stopcondition(void *Apointer, double *x, double *b)
{
  MX *A;
  A = Apointer;

  if ( B != 0.0 ) {
    double norm;
    norm = psc98(A,x,b);
    if ( defop("-v") ) fprintf(stderr,"%e\n",norm);
    if ( norm < B ) return 1;
    else return 0;
  }
  return 1;
}


int estiva_psc98condition(double *x, double *b)
{
  return stopcondition(A,x,b);
}
