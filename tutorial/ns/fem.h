#ifndef _ESTIVA_FEM_H_
#define _ESTIVA_FEM_H_

#include <stdio.h>
#include "estiva/mesh.h"
#include "estiva/que.h"

#define femdelta(S,Z,N)  estiva_femdelta(&S,Z,N)
#define dimp2(N)         estiva_dimp2(N)
#define rectmesh(Z,N) estiva_rectmesh(&Z,&N)
#define pltp2(x,Z,N)  estiva_pltp2(x, Z, N)
#define forgammap2(i,label,Z,N)						\
  for( estiva_forqinit(estiva_forgammap2_init(Z,N,label),		\
		       (void*)&estiva_forgammap2_p);			\
       estiva_forgammap2_flag = estiva_forq((void*)&estiva_forgammap2_p), \
	 i = *estiva_forgammap2_p,					\
	 estiva_forgammap2_flag;					\
       )

void estiva_femdelta(double **Sp, xyc *Z, nde *N);
long estiva_dimp2(nde *N);
void estiva_rectmesh(xyc **Zp, nde **Np);
void estiva_pltp2(double *x, xyc * Z, nde *N);
que *estiva_forgammap2_init(xyc *Z, nde *N, char *label);


extern long *estiva_forgammap2_p, estiva_forgammap2_flag;

#endif
