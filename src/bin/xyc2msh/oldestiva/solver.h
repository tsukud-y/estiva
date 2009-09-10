#ifndef _SOLVER_H_
#define _SOLVER_H_
extern int solver_gauss(void *A, void *B);
extern int solver_inv(void *A);

#define gauss(A,B) solver_gauss(A,B)
#define inv(A)     solver_inv(A)
#endif
