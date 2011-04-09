#ifndef _ESTIVA_SOLVER_H_
#define _ESTIVA_SOLVER_H_

void   *estiva_distmx(void *A);
double *estiva_matvecmpi(void *A,double alpha,double *p,double beta,double *q);
long    estiva_sendcommand(long command);
int     estiva_solver(void* A, double *x, double* b);
int     estiva_MPI_np;

#define distmx(Apointer)                estiva_distmx(Apointer)
#define matvecmpi(A, alpha, p, beta, q) estiva_matvecmpi(A,alpha,p,beta,q)
#define sendcommand(command)            estiva_sendcommand(command)
#define solver(A, x, b) estiva_solver(A, x, b)

#endif
