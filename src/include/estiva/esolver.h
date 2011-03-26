#ifndef _ESTIVA_ESOLVER_H_
#define _ESTIVA_ESOLVER_H_

double  estiva_esolver(double **A, double *x);

#define esolver(A, x) estiva_esolver(A, x)

#endif
