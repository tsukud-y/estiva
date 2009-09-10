#ifndef _ESTIVA_SOLVER_H_
#define _ESTIVA_SOLVER_H_

#define solver(pA, x, b) estiva_solver(pA, x, b)

extern int estiva_solver(void* pA, double *x, double* b);

#endif
