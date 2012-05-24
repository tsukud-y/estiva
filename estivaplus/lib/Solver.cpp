#include "estivaplus.h"
#include "viennacl.h"

vector gpusolver( matrix A, vector b)
{
  int n = A.size();
  vector x(n);
  gpumatrix Agpu(n,1); gpuvector bgpu(n), xgpu(n);

  copy(A, Agpu); copy(b.begin(), b.end(), bgpu.begin());

  ILU     vcl_ilut(Agpu, ilut_tag(8,1e-3));
  //Scaling vcl_row_scaling(Agpu, row_scaling_tag(2));
  //Jacobi  vcl_jacobi(Agpu,jacobi_tag());


  bicgstab_tag  custom_cg(5e-7,1000);
  xgpu = solve(Agpu, bgpu, custom_cg, vcl_ilut);
  //xgpu = solve(Agpu, bgpu, custom_cg, vcl_jacobi);
  //bicgstab_tag  custom_bicgstab(5e-7,1000000);
  //xgpu = solve(Agpu, bgpu, custom_bicgstab, vcl_jacobi);

  copy(xgpu.begin(), xgpu.end(), x.begin());
  return x;
}

extern vector gpusolver(matrix A, vector b);

void Solver(Matrix &Am, Vector &xv, Vector &bv)
{
  xv = gpusolver(Am,bv);
  return;
  Matrix B(2); Vector x(2), b(2);
  
  B[0][0] = 1.0; B[0][1] = 0.5; 
  B[1][0] = 0.0; B[1][1] = 1.0; 


  b[0] = 1.0;
  b[1] = 1.0;


  x = gpusolver(B,b);

  int i;

  i = 0; printf("x[%d] = %f\n",i,x[i]);
  i = 1; printf("x[%d] = %f\n",i,x[i]);

}

void Solverorg(Matrix &Am, Vector &xv, Vector &bv)
{
  static MX *A;
  static double *b, *x;
  long NUM = Am.capacity();

  ary1(b,NUM+1);
  ary1(x,NUM+1);
  Matrix2mx(Am,&A);
  Vector2ary(bv,b);
  estiva_qmrsolver(A,x,b);
  xv = ary2Vector(x);
}




#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>

#define VIENNACL_HAVE_UBLAS 1

#define SOLVER_ITERS 2500
//#define SCALAR float
#define SCALAR double

//#define SOLVER_TOLERANCE 1e-5
#define SOLVER_TOLERANCE 1e-9

#include "viennacl/vector.hpp"
#include "viennacl/coordinate_matrix.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/io/matrix_market.hpp"
#include "viennacl/linalg/norm_2.hpp"

#include "viennacl/linalg/amg.hpp"

#include <iostream>
#include <vector>
#include <ctime>
#include "vector-io.hpp"


extern void run_amg(viennacl::linalg::bicgstab_tag & cg_solver,
             viennacl::vector<double> & vcl_vec,
             viennacl::vector<double> & vcl_result,
             viennacl::compressed_matrix<double> & vcl_compressed_matrix,
             viennacl::linalg::amg_tag & amg_tag);


int main2()
{
  
  typedef double ScalarType; 

  boost::numeric::ublas::vector<ScalarType> ublas_vec, ublas_result;
  boost::numeric::ublas::compressed_matrix<ScalarType> ublas_matrix;
  
  viennacl::linalg::bicgstab_tag cg_solver;
  viennacl::linalg::amg_tag amg_tag;
  viennacl::linalg::amg_precond<boost::numeric::ublas::compressed_matrix<ScalarType> > ublas_amg;
    
  viennacl::io::read_matrix_market_file(ublas_matrix, "../examples/testdata/mat65k.mtx");
  
  readVectorFromFile("../examples/testdata/rhs65025.txt", ublas_vec);
  readVectorFromFile("../examples/testdata/result65025.txt", ublas_result);
  
  viennacl::vector<ScalarType> vcl_vec(ublas_vec.size());
  viennacl::vector<ScalarType> vcl_result(ublas_vec.size());
  viennacl::compressed_matrix<ScalarType> vcl_compressed_matrix(ublas_vec.size(), ublas_vec.size());

  // Copy to GPU
  viennacl::copy(ublas_matrix, vcl_compressed_matrix);
  viennacl::copy(ublas_vec, vcl_vec);
  viennacl::copy(ublas_result, vcl_result);

  amg_tag = viennacl::linalg::amg_tag(VIENNACL_AMG_COARSE_RS,
                                      VIENNACL_AMG_INTERPOL_DIRECT,
                                      0.25,
                                      0.2, 
                                      0.67,
                                      3,   
                                      3,   
                                      0);  
  run_amg (cg_solver,  vcl_vec, vcl_result, vcl_compressed_matrix, amg_tag);
  
  return 0;
}

vector gpusolveramg( matrix A, vector b)
{
  int n = A.size();
  vector x(n);
  gpumatrix Agpu(n,1); gpuvector bgpu(n), xgpu(n);

  copy(A, Agpu); copy(b.begin(), b.end(), bgpu.begin());


  viennacl::linalg::bicgstab_tag bicgstab_solver(5e-7,1000);
  viennacl::linalg::amg_tag amg_tag;

  switch(atoi(getop("-case"))){
  case 1: amg_tag = viennacl::linalg::amg_tag(VIENNACL_AMG_COARSE_RS, VIENNACL_AMG_INTERPOL_DIRECT, 0.25, 0.2, 0.67, 3,3, 0); break;
  case 2: amg_tag = viennacl::linalg::amg_tag(VIENNACL_AMG_COARSE_RS, VIENNACL_AMG_INTERPOL_CLASSIC, 0.25, 0.2, 0.67, 3, 3, 0); break;
  case 3: amg_tag = viennacl::linalg::amg_tag(VIENNACL_AMG_COARSE_ONEPASS, VIENNACL_AMG_INTERPOL_DIRECT,0.25, 0.2, 0.67, 3, 3, 0); break;
  case 4: amg_tag = viennacl::linalg::amg_tag(VIENNACL_AMG_COARSE_RS0, VIENNACL_AMG_INTERPOL_DIRECT, 0.25, 0.2, 0.67, 3, 3, 0); break;
  case 5: amg_tag = viennacl::linalg::amg_tag(VIENNACL_AMG_COARSE_RS3, VIENNACL_AMG_INTERPOL_DIRECT, 0.25, 0.2, 0.67, 3, 3, 0); break;
  case 6: amg_tag = viennacl::linalg::amg_tag(VIENNACL_AMG_COARSE_AG, VIENNACL_AMG_INTERPOL_AG, 0.08, 0, 0.67, 3, 3, 0); break;
  case 7: amg_tag = viennacl::linalg::amg_tag(VIENNACL_AMG_COARSE_AG, VIENNACL_AMG_INTERPOL_SA, 0.08, 0.67, 0.67, 3, 3, 0); break;
  }


  run_amg (bicgstab_solver, bgpu, xgpu, Agpu, amg_tag);

  copy(xgpu.begin(), xgpu.end(), x.begin());
  return x;
}

