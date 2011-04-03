#include <stdio.h>
#include <stdlib.h>
#include <estiva/ary.h>
#include <estiva/mx.h>

#define rmx(A,i,j)            estiva_rmx(A,i,j)
#define wmx(A,i,j,a)          estiva_wmx(A,i,j,a)

void estiva_initmx(MX **Ap, long i, long j){
  MX *T;
  if(*Ap == NULL) ary1(*Ap,1); //*Ap = calloc(1,sizeof(MX));
  T = *Ap;
  T->m = i-1;
  T->n = j-1;
  T->I = 1;
  T->J = 1;
  ary2(T->A, T->m, T->n);
  ary2(T->IA,T->m, T->n);
}

static double estiva_rmx(MX *T, long i, long j){
  register double **A;
  register long   J,  *IAi, **IA, n;
  A  = T->A;
  IA = T->IA;
  n  = T->n;
  IAi = IA[i-1];
  for(J=0; J<n; J++){
    if(IAi[J]==0) return       0.0;
    if(IAi[J]==j) return A[i-1][J];
  } 
  return 0.0;
}

static void estiva_wmx(MX *T, long i, long j, double a){
  register double **A;
  register long  J, **IA, *IAi, IAiJ, n;
  A  = T->A;
  IA = T->IA;
  n  = T->n;
  if ( a == rmx(T,i,j) ) return;
  IAi = IA[i-1];
  for(J=0; J<n; J++){
    IAiJ = IAi[J];
    if (IAiJ == 0 || IAiJ ==j){
      A[i-1][J] = a;
      if (IAiJ==0) IAi[J] = j; 
      return;
    }
  } 
  fprintf(stderr,"initmx(m,n)'s n is too short\n");
  abort();
}

double *estiva_mx(MX *T, long i, long j){

  wmx(T,T->I,T->J,T->a);
  T->I = i; 
  T->J = j;
  T->a = rmx(T,i,j);
  return &(T->a);
}
