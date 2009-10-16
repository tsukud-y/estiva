#ifndef _ESTIVA_PRECOND_H_
#define _ESTIVA_PRECOND_H_

#define precondnone(n, x, b) estiva_precondnone(n, x, b)
#define precondscaling(x, D, b) estiva_precondscaling(x, D, b)

extern void estiva_precondnone(long n, double *x, double* b);
extern void estiva_precondscaling(double *x, double *D, double* b);

#endif
