#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "estiva/op.h"
#include "estiva/mx.h"
#include "estiva/std.h"
#include "estiva/ary.h"
#include "estiva/vec.h"

static int estiva_MPI_np;
static long estiva_mpicommand(long command);
static MX *estiva_A;
static int estiva_efenceflag;

void estiva_efence(int flag)
{
  estiva_efenceflag = flag;
}


void estiva_distributemx(void *Apointer)
{
  MX *A;
  estiva_A = A = Apointer;

  if (defop("--mpi")) estiva_efence(1);
  if (defop("-mpi")) estiva_efence(0);

  if ( estiva_efenceflag ) return;

  mx(A,1,1) = mx(A,1,1);
  estiva_mpicommand(1);
  MPI_Bcast(&A->n, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A->w, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A->A[0][0], A->n*A->w, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A->IA[0][0], A->n*A->w, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A->I, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A->J, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A->a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

static void estiva_recievemx(MX **Apointer, int p)
{
  MX *A;
  long i, j, n, w;

  MPI_Bcast(&n, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&w, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  initmx(*Apointer,n+1,w+1);
  A = *Apointer;

  MPI_Bcast(&A->A[0][0], A->n*A->w, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A->IA[0][0], A->n*A->w, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A->I, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A->J, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A->a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  forall (0, j, A->w-1) forall (0, i, A->n-1) { 
    if ((p-1)*A->n/(estiva_MPI_np-1)<=i && i<p*A->n/(estiva_MPI_np-1)) {
    } else {
      mx(A,i+1,A->IA[i][j]) = 0.0;
    }
  }

#if 0
  if ( p == 2) {
    printf("%d:A(1,1)=%f A(1,2)=%f A(1,3)=%f \n",p,mx(A,1,1),mx(A,1,2),mx(A,1,3));
    printf("%d:A(2,1)=%f A(2,2)=%f A(2,3)=%f \n",p,mx(A,2,1),mx(A,2,2),mx(A,2,3));
    printf("%d:A(3,1)=%f A(3,2)=%f A(3,3)=%f \n",p,mx(A,3,1),mx(A,3,2),mx(A,3,3));
  }
#endif


}

static void slaveserver(int p)
{
  static MX *A;
  static double  alpha, *x, beta, *y;
  long command;

  while(1){
    MPI_Bcast(&command, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    switch(command) {
    case 1:
      estiva_recievemx((void*)&A,p);
      break;
    case 2:
      ary1(x,getveclength());
      ary1(y,getveclength());
      MPI_Bcast(&alpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(x,getveclength(),MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&beta, 1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(y,getveclength(),MPI_DOUBLE,0,MPI_COMM_WORLD);
      matvecvec(A,alpha,x,beta,y);
      MPI_Reduce(y,NULL,getveclength(),MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
      break;
    default:
      MPI_Finalize();
      if (defop("-v")) fprintf(stderr,"Finazlie %d \n",p);
      exit(0);
    }
  }
}


static void estiva_atexitmpimaster(void)
{
  estiva_mpicommand(65535);
  MPI_Finalize();
  if (defop("-v"))fprintf(stderr,"Finalize 0 \n");
}


static int estiva_initmpi(void)
{
  static int init;
  int p;
  
  if ( init == 0 ){
    init = 1;
    MPI_Init(estiva_pargc, estiva_pargv);
    MPI_Comm_size(MPI_COMM_WORLD,&estiva_MPI_np);
    MPI_Comm_rank(MPI_COMM_WORLD,&p);
    if ( p == 0 ) atexit(estiva_atexitmpimaster);
    else  slaveserver(p);
  }
  return estiva_MPI_np;
}

static long estiva_mpicommand(long command)
{
  estiva_initmpi();
  MPI_Bcast(&command, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  return command;
}

static void estiva_matvecmpi(alpha,x,beta,y)
     double alpha, *x, beta, *y;
{
  static double *tmp;
  ary1(tmp,getveclength());
  estiva_mpicommand(2);
  MPI_Bcast(&alpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(x,getveclength(),MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&beta, 1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(y,getveclength(),MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Reduce(tmp,y,getveclength(),MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
}

int estiva_matvecmpi2(alpha,x,beta,y)
     double *alpha, *x, *beta, *y;
{
  if ( estiva_efenceflag || estiva_MPI_np <= 1)
    matvecvec(estiva_A,*alpha,x,*beta,y);
  else
    estiva_matvecmpi(*alpha,x,*beta,y);
  return 0;
}


