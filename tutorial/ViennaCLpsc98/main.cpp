#include <stdio.h>

/*#include "viennacl.h"                                 */
#ifndef _VIENNACL_H_
#define _VIENNACL_H_
#include <iostream>
#include <viennacl/scalar.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/compressed_matrix_operations.hpp>
#include <viennacl/linalg/row_scaling.hpp>
#include <viennacl/linalg/jacobi_precond.hpp>

using viennacl::linalg::ilut_precond;
using viennacl::linalg::jacobi_precond;
using viennacl::linalg::row_scaling;
using viennacl::compressed_matrix;
using viennacl::linalg::cg_tag;
using viennacl::linalg::solve;
using viennacl::linalg::ilut_tag;
using viennacl::linalg::jacobi_tag;
using viennacl::linalg::row_scaling_tag;

typedef compressed_matrix<double>  SparseMatrix;
typedef std::vector< std::map< unsigned int, double> > matrix;
typedef std::vector<double> vector;
typedef viennacl::compressed_matrix<double> gpumatrix;
typedef viennacl::vector<double> gpuvector;
typedef ilut_precond< SparseMatrix >  ILU;
typedef row_scaling< SparseMatrix > Scaling;
typedef jacobi_precond< SparseMatrix > Jacobi;

extern "C" {
  int genmat(int,int*,double*,double*);
  int chkval(FILE*,int,double*);
}

vector gpusolver( matrix A, vector b)
{
  int n = A.size();
  vector x(n);
  gpumatrix Agpu(n,1); gpuvector bgpu(n), xgpu(n);

  copy(A, Agpu); copy(b.begin(), b.end(), bgpu.begin());

  ILU     vcl_ilut(Agpu, ilut_tag(8,1e-3));
  Scaling vcl_row_scaling(Agpu, row_scaling_tag(2));
  Jacobi  vcl_jacobi(Agpu,jacobi_tag());
  cg_tag  custom_cg(1e-7,1000000);

  xgpu = solve(Agpu, bgpu, custom_cg, vcl_jacobi);

  copy(xgpu.begin(), xgpu.end(), x.begin());
  return x;
}
#endif
/***********************************************************/


int main(int argc, char **argv){
  vector AA(10);
  double B;
  int    JA[10], i, j, n, w;

  genmat(-1,&JA[0],&AA[0],&B);
  n = JA[0]; w = JA[2];

  matrix A(n); vector x(n), b(n);

  for ( i=1; i<=n; i++) {
    for (j=0;j<=w-1;j++) {
      JA[j] =  -1;
      AA[j] = 0.0;
    }
    genmat(i,&JA[0],&AA[0],&B);
    for ( j=0; j<w; j++) if (JA[j] != -1){
	if( JA[j] <=0 || n < JA[j] ) {
	  ;
	} else {
	  A[i-1][JA[j]-1] = AA[j];
	  b[i-1] = B;
	}
      }
  }

  x = gpusolver(A,b);
  
  chkval(stdout,n,&x[0]);
  return 0;
}
