#ifndef _ESTIVA_FEM_H_
#define _ESTIVA_FEM_H_

#include <stdio.h>
#include "estiva/mesh.h"
#include "estiva/que.h"

#define femdelta(S,Z,N)  estiva_femdelta(&S,Z,N)
void estiva_femdelta(double **Sp, xyc *Z, nde *N);
void setZNS(xyc *Z, nde *N, double *S);

#endif
