#ifndef _ESTIVA_FEM_H_
#define _ESTIVA_FEM_H_

#include <stdio.h>
#include "estiva/mesh.h"
#include "estiva/que.h"

#define femdelta(S,Z,N)  estiva_femdelta(&S,Z,N)
#define forgammap2(i,label,Z,N)						\
  for( estiva_forqinit(estiva_forgammap2_init(Z,N,label),		\
		       (void*)&estiva_forgammap2_p);			\
       estiva_forgammap2_flag = estiva_forq((void*)&estiva_forgammap2_p), \
	 i = *estiva_forgammap2_p,					\
	 estiva_forgammap2_flag;					\
       )

void estiva_femdelta(double **Sp, xyc *Z, nde *N);
que *estiva_forgammap2_init(xyc *Z, nde *N, char *label);
void setZNS(xyc *Z, nde *N, double *S);


extern long *estiva_forgammap2_p, estiva_forgammap2_flag;

#endif
