#include <stdio.h>
#include <stdlib.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/op.h"
#include "estiva/vec.h"
#include "estiva/std.h"
#include "/opt/mpich2/include/mpi.h"
#include "estiva/solver.h"

static long dim1vec;
static int np;

void *estiva_distmx(void *Apointer)
{
  static MX *A;
  static double *AA;
  static long   *IA;
  long n, i, j, p, w;


  MPI_Comm_size(MPI_COMM_WORLD,&np);
  np--;

  A = Apointer;
  mx(A,1,1) = mx(A,1,1);
  n = A->m;
  w = A->n;
  ary1(AA,n+1);
  ary1(IA,n+1);
  dim1vec = n;

  for (p = 1; p<=np; p++) {
    MPI_Send(&w,1,MPI_LONG,p, 999, MPI_COMM_WORLD);
    MPI_Send(&n,1,MPI_LONG,p  ,1000,  MPI_COMM_WORLD);
  }

  forall (1,p,np) forall(0, j, w-1) {
    forall(0, i, n-1) { 
      if ( (p-1)*n/np  <= i && i < p*n/np ) {
	AA[i] = A->A[i][j];
	IA[i] = A->IA[i][j];
      } else {
	AA[i] = 0.0;
	IA[i] = 0;
      }
    }
    MPI_Send(AA,n,MPI_DOUBLE,p,1001+j,MPI_COMM_WORLD);
    MPI_Send(IA,n,MPI_LONG,  p,2001+j,MPI_COMM_WORLD);
  }
  return A;
}

static void *recvmx(void **Appointer)
{
  static MX     *A, **App;
  static double *AA;
  static long   *IA;
  MPI_Status status;
  long  i, j, n, w;

  MPI_Recv(&w, 1, MPI_LONG,   0,  999, MPI_COMM_WORLD, &status);
  MPI_Recv(&n, 1, MPI_LONG,   0, 1000, MPI_COMM_WORLD, &status);
  dim1vec = n;
  App = (MX **)Appointer;
  initmx(*App,n+1,w+1);
  A = *App;
  clearmx(A);
  ary1(AA,n);
  ary1(IA,n);
  forall(0,j,w-1){
    MPI_Recv(AA, n, MPI_DOUBLE, 0, 1001+j, MPI_COMM_WORLD, &status);
    MPI_Recv(IA, n, MPI_LONG,   0, 2001+j, MPI_COMM_WORLD, &status);
    forall (0,i,n-1) mx(A,i+1,IA[i]) = AA[i];
  }
  setveclength(n);
  return A;
}

static void *distvec(double *x)
{
  MPI_Bcast(x,dim1vec,MPI_DOUBLE,0,MPI_COMM_WORLD);
  return x;
}

static void *recvvec(double *x)
{
  MPI_Bcast(x,dim1vec,MPI_DOUBLE,0,MPI_COMM_WORLD);
  return x;
}
 
static void *returnvec(double *x)
{
  MPI_Reduce(x,NULL,dim1vec,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
  return x;
}

static void *collectvec(double *y)
{
  static double *x;
  long i;

  ary1(x,dim1vec);
  forall (0,i,dim1vec-1) x[i] = 0.0;
  forall (0,i,dim1vec-1) y[i] = 0.0;  
  MPI_Reduce(x,y,dim1vec,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
  return y;
}

static void slave(int p)
{
  static MX *A;
  static double *x,*y;
  long command;
  
  while ( 1 ) {
    MPI_Bcast(&command, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    if (command == 1) { recvmx((void*)&A);}
    if (command == 2) {
      ary1(x,A->m);
      ary1(y,A->m);
      setveclength(A->m);
      recvvec(x);
      matvecvec(A,1.0,x,0.0,y);
      returnvec(y);
    }
    if (command == 999) {MPI_Finalize(); exit(0);}
  }
}

static int mpisolverinit(void)
{
  static int init;
  int p;

  if ( init == 0 ){
    init = 1;
    MPI_Init(estiva_pargc, estiva_pargv);
    MPI_Comm_rank(MPI_COMM_WORLD,&p);
    if ( p != 0 ) {slave(p); exit(0);}
  }
  return init;
}

long estiva_sendcommand(long command)
{
  mpisolverinit();
  MPI_Bcast(&command, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  if ( command == 999) { MPI_Finalize(); exit(0);}
  return command;
}

double *estiva_matvecmpi(void *Apointer, double alpha, double *p, 
		  double beta, double *q)
{
  MX *A;
  long i;
  A = Apointer;

  forall (0,i,A->m-1) q[i] = alpha*p[i];
  sendcommand(2);
  distvec(q);
  collectvec(q);

  return q;
}
