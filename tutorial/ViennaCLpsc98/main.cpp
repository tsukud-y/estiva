//#include <stdio.h>
//#include "estiva/ary.h"
//#include "estiva/mx.h"
//#include "estiva/op.h"
//#include "estiva/solver.h"
//#include "estiva/std.h"
//#include "estiva/vec.h"
//#include "estiva/eblas.h"
#include <iostream>
#include <viennacl/scalar.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/compressed_matrix_operations.hpp>
#include <viennacl/linalg/row_scaling.hpp>
#include <viennacl/linalg/jacobi_precond.hpp>



extern "C" {
int genmat(int,int*,double*,double*);
int chkval(FILE*,int,double*);
}

int main(int argc, char **argv){
  //static MX *A;
  static int JA[10];
  static double AA[10], b[2000001], x[2000001];
  double B;
  int    i, j, n, w;

  //initop(argc,argv);


  genmat(-1,JA,AA,&B);
  n = JA[0]; w = JA[2];
  std::vector< std::map< unsigned int, double> > cpu_sparse_matrix(n);
  std::vector<double> stl_vec(n);


  //initmx(A,n+1,w+1);  ary1(x,n+1);  ary1(b,n+1);

  for ( i=1; i<=n; i++) {
    for (j=0;j<=w-1;j++) {
      JA[j] =  -1;
      AA[j] = 0.0;
    }
    genmat(i,JA,AA,&B);
    for ( j=0; j<w; j++) if (JA[j] != -1){
	if( JA[j] <=0 || n < JA[j] ) {
	  ;
	} else {
	  //mx(A,i,JA[j]) = AA[j]; 
	  cpu_sparse_matrix[i-1][JA[j]-1] = AA[j];
	  stl_vec[i-1] = b[i] = B;
	}
      }
  }
  viennacl::compressed_matrix<double> vcl_sparse_matrix(n,1);
  viennacl::vector<double> vcl_vec(n);
  viennacl::vector<double> vcl_result(n);
  copy(cpu_sparse_matrix, vcl_sparse_matrix);
  copy(stl_vec.begin(), stl_vec.end(), vcl_vec.begin());


  viennacl::linalg::cg_tag custom_cg(1e-7,100000);



  using viennacl::compressed_matrix;

  typedef compressed_matrix<double>  SparseMatrix;
  
  //ilut_precond<SparseMatrix> vcl_ilut(vcl_sparse_matrix, viennacl::linalg::ilut_tag());


  using viennacl::linalg::ilut_precond;
  using viennacl::linalg::jacobi_precond;
  using viennacl::linalg::row_scaling;

#if 0
  ilut_precond< SparseMatrix > vcl_ilut(vcl_sparse_matrix,
					viennacl::linalg::ilut_tag(8,1e-3));


  row_scaling< SparseMatrix > vcl_row_scaling(vcl_sparse_matrix,
					      viennacl::linalg::row_scaling_tag(2));
#endif

  jacobi_precond< SparseMatrix > vcl_jacobi(vcl_sparse_matrix,
					    viennacl::linalg::jacobi_tag());


  vcl_result = viennacl::linalg::solve(vcl_sparse_matrix,vcl_vec,
				       custom_cg, vcl_jacobi);


  copy(vcl_result.begin(), vcl_result.end(), stl_vec.begin());

  for ( i = 0; i<n; i++) x[i] = stl_vec[i];

  chkval(stdout,n,&x[0]);


  return 0;
}
