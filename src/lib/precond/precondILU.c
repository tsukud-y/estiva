#include <stdio.h>
#include <stdlib.h>
#include <estiva/ary.h>
#include <estiva/solver.h>
#include <estiva/mx.h>
void fprintCRS(FILE *fp, CRS crs);

void mx2CRS(MX *A, CRS *crs)
{
  long i, j, k, n=0, *diag_ptr, *row_ptr, *col_ind,found, FALSE=0, TRUE=1;
  double *val, *pivots, element=0.0;

  mx(A,1,1) = mx(A,1,1);

  for ( i = 1;  i <= A->m;  i++  )  
    for (  j = 1;  j <= A->m;  j++  )
      if ( mx(A,i,j) != 0.0 || i == j) n++;

  crs->m = A->m;
  crs->n = n;
  ary1(crs->val, n+1);
  ary1(crs->col_ind, n+1);
  ary1(crs->row_ptr, A->m + 2);
  ary1(crs->diag_ptr, A->m + 1);
  ary1(crs->pivots, A->m +1);

  k = 1;
  for ( i = 1;  i <= A->m;  i++  )  {
    crs->row_ptr[i] = k;
    for (  j = 1;  j <= A->m;  j++  )  {
      if ( mx(A,i,j) != 0.0 || i == j) {
	crs->val[k] = mx(A,i,j);
	crs->col_ind[k] = j;
	if (i == j) crs->diag_ptr[i] = k;
	k++;
      }
    }
  }
  crs->row_ptr[i] = k;



  /* Incomplete LU Decomposition */
  n        = A->m;
  val      = crs->val;
  diag_ptr = crs->diag_ptr;
  col_ind  = crs->col_ind;
  row_ptr  = crs->row_ptr;
  pivots   = crs->pivots;

  for ( i = 1;  i <= n;  i++ )  {
    pivots[i] = val[diag_ptr[i]];
  }

  for  ( i = 1;  i <=n;  i++ )  {
    pivots[i] = 1.0/pivots[i];
    for  ( j = diag_ptr[i]+1;  j <= row_ptr[i+1]-1;  j++)  {
      found = FALSE;
      for  (  k = row_ptr[col_ind[j]];  k <= diag_ptr[col_ind[j]]-1;  k++)  {
	if  (  col_ind[k] == i  )  {
	  found = TRUE;
	  element = val[k];
	}  // endif
      }  // end;
      
      if  (  found == TRUE  )
	val[diag_ptr[col_ind[j]]] = val[diag_ptr[col_ind[j]]]
	  - element * pivots[i] * val[j];
    }
  }
}


void solveILU(CRS crs, double *y, double *x)
{
  long i, j, n, *row_ptr, *diag_ptr, *col_ind;
  double sum, *val, *pivots;
  static double *z;

  n        = dim1(crs.pivots);
  row_ptr  = crs.row_ptr;
  diag_ptr = crs.diag_ptr;
  val      = crs.val;
  col_ind  = crs.col_ind;
  pivots   = crs.pivots;

  ary1(z, n+1);
  x--, y--;

  for  (  i = 1;  i <= n;  i++  )  { y[i] = z[i] = 0.0; }

  for  (  i = 1;  i <= n;  i++  )  {
    sum = 0.0;
    for  (  j = row_ptr[i];  j <= diag_ptr[i]-1;  j++  )  {
      sum = sum + val[j] * z[col_ind[j]];
    }
    z[i] = pivots[i] * ( x[i] - sum );
  }

  for  (  i = n;  i >= 1;  i--  )  {
    sum = 0.0;
    for  (  j = diag_ptr[i];  j <= row_ptr[i+1]-1;  j++  )  {
      sum = sum + val[j] * y[col_ind[j]];
      y[i] = z[i] - pivots[i] * sum;
    }
  }
}


void solveILUT(CRS crs, double *y, double *x)
{
  long i, j, n, *row_ptr, *diag_ptr, *col_ind;
  double *val, *pivots, tmp;
  static double *z, *x_tmp;

  n        = dim1(crs.pivots);
  row_ptr  = crs.row_ptr;
  diag_ptr = crs.diag_ptr;
  val      = crs.val;
  col_ind  = crs.col_ind;
  pivots   = crs.pivots;

  ary1(z, n+1);
  ary1(x_tmp, n+1);


  for  (  i = 1;  i <= n;  i++  )  {
    z[i] = 0.0;
    x_tmp[i] = x[i];
  }

  for  (  i = 1;  i <= n;  i++  )  {
    z[i] = x_tmp[i];
    tmp = pivots[i] * z[i];
    for  (  j = diag_ptr[i]+1;  j <= row_ptr[i+1]-1;  j++  )  {
      x_tmp[col_ind[j]] = x_tmp[col_ind[j]] - tmp * val[j];
    }
  }

  for  (  i = n;  i >= 1;  i--  )  {
    y[i] = pivots[i] * z[i];
    for  (  j = row_ptr[i];  j <= diag_ptr[i]-1;  j++  )  {
      z[col_ind[j]] = z[col_ind[j]] - val[j] * y[i];
    }
  }
} 


void fprintCRS(FILE *fp, CRS crs)
{
  long k, m, n;
  n = crs.n;
  m = crs.m;

  fprintf(fp,"val ");
  for ( k = 1;  k<=n;  k++) {
    fprintf(fp,"%.1f ",crs.val[k]);
  }
  fprintf(fp,"\n");

  fprintf(fp,"col_ind ");
  for ( k = 1;  k<=n;  k++) {
    fprintf(fp,"%ld ",crs.col_ind[k]);
  }
  fprintf(fp,"\n");

  fprintf(fp,"row_ptr ");
  for ( k = 1;  k<= crs.m+1;  k++) {
    fprintf(fp,"%ld ",crs.row_ptr[k]);
  }
  fprintf(fp,"\n");

  fprintf(fp,"diag_ptr ");
  for ( k = 1;  k<= crs.m;  k++) {
    fprintf(fp,"%ld ",crs.diag_ptr[k]);
  }
  fprintf(fp,"\n");

  fprintf(fp,"pivots ");
  for ( k = 1;  k<= crs.m;  k++) {
    fprintf(fp,"%.1f ",crs.pivots[k]);
  }
  fprintf(fp,"\n");
}



void estiva_precondILU(CRS *dummypivots, MX *A, double *x, double *b)
{
  solveILU(*dummypivots,x,b);
}

void estiva_ILU(CRS *dummypivots, MX *A)
{
  long i;
  
  for  (  i = 1;  i < A->m;  i++ ) if ( mx(A,i,i) == 0.0 ) mx(A,i,i) = 1.0;
  
  mx2CRS(A,dummypivots);

  //if (init == 1) {fprintCRS(fopen("CRST.txt","w"), *dummypivots); exit(0);}
  //init = 1;
}
