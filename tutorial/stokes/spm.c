#include <stdio.h>
#include <stdlib.h>
#include "ary.h"
#include "spm.h"

typedef struct{
  long   j;
  double x;
  void*  p;
} row;

typedef struct{
  long   i, j;
  double x;
} cache;

static row *Ai;

static double spm_read(void *i, long j)
{
  for(Ai=i;Ai;Ai=Ai->p) if(Ai->j == j) return Ai->x;
  return 0.0;
}

static void spm_write(spm* A, long i, long j, double x)
{
  if(A[i] == NULL) A[i] = calloc(sizeof(row),1);
  
  for(Ai=(void*)A[i];Ai->p;Ai=Ai->p) if(Ai->j == j){ Ai->x = x; return; }
  if(x == 0.0) return;
  Ai->j = j; 
  Ai->x = x; 
  Ai->p = calloc(sizeof(row),1);
  return;
}

double* spm_double(spm* A, long i, long j)
{
  static cache *c;

  if(A[0] == NULL) A[0] = calloc(sizeof(cache),1);
  c  = A[0];
  spm_write(A,c->i,c->j,c->x);

  c->i = i;
  c->j = j; 
  c->x = spm_read(A[i],j);

  return &(c->x);
}

static int comp(a, b)
     long *a, *b;
{
  if( *a < *b ) return  1;
  if( *a > *b ) return -1;
  return 0;
}

long* nonzero_col(spm *p)
{
  static long n;
  static long *J;
  
  for(n=1,Ai=(void*)p;Ai->p;Ai=Ai->p) n++;
  ary1(J,n);
  for(n=1,Ai=(void*)p;Ai->p;Ai=Ai->p) J[n++] = Ai->j;

  qsort(&J[1],dim1(J),dim0(J),comp);
  return J;
}
