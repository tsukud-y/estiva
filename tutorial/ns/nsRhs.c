#include "ns.h"
#include "fem.h"


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include "estiva/mesh.h"
#include "estiva/ary.h"
#include "estiva/op.h"
#include "estiva/mx.h"
#include "estiva/foreach.h"
#include "estiva/que.h"
#include "estiva/solver.h"


void estiva_nsRhs(double *b, xyc *Z, nde *N, MX *M, double *x, double t)
{
  long   i, j, NUM, m, n;

  m = dimp2(N);
  n = dim1(Z);

  NUM = m*2+n;

  for( i = 1; i <= NUM; i++ ) b[i] = 0.0;

  for( i = 1; i <= m; i++ ) for ( j = 1; j <= m; j++ ) {
      b[  i] += mx(M,i,j)*x[  j];
      b[m+i] += mx(M,i,j)*x[m+j];
    }


  /*                                                                            
  static int init = 0;                                                          
  static MX *M2;                                                                
                                                                                
  if ( !init ) {                                                                
    long i, j, m, n;                                                            
    init = 1;                                                                   
    m = dimp2(N);                                                               
    n = dim1(Z);                                                                
                                                                                
    initmx(M2, 2*m+n, 50);                                                      
                                                                                
    for ( i = 1; i <= m; i++ ) for ( j = 1; j <=m; j++ ) {                      
        mx(M2,  i,  j) = mx(M,i,j);                                             
        mx(M2,m+i,m+j) = mx(M,i,j);                                             
      }                                                                         
  }                                                                             
  mulmx(b,M2,x);                                                                
  */
}
