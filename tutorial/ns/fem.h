#ifndef _ESTIVA_FEM_H_
#define _ESTIVA_FEM_H_

#include <stdio.h>
#include "estiva/mesh.h"

#define femdelta(S,Z,N)  estiva_femdelta(&S,Z,N)
#define dimp2(N)         estiva_dimp2(N)
#define rectmesh(Z,N) estiva_rectmesh(&Z,&N)


void estiva_femdelta(double **Sp, xyc *Z, nde *N);
long estiva_dimp2(nde *N);
void estiva_rectmesh(xyc **Zp, nde **Np);


#endif
