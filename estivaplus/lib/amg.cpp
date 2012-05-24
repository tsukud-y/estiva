/* =========================================================================
   Copyright (c) 2010-2011, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at
               
   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

#ifndef NDEBUG     //without NDEBUG the performance of sparse ublas matrices is poor.
 #define NDEBUG
#endif

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


template <typename MatrixType, typename VectorType, typename SolverTag, typename PrecondTag>
void run_solver(MatrixType const & matrix, VectorType const & rhs, VectorType const & ref_result, SolverTag const & solver, PrecondTag const & precond)
{
  VectorType result(rhs);
  VectorType residual(rhs);
  
  result = viennacl::linalg::solve(matrix, rhs, solver, precond);
  residual -= viennacl::linalg::prod(matrix, result);
  std::cout << "  > Relative residual: " << viennacl::linalg::norm_2(residual) / viennacl::linalg::norm_2(rhs) << std::endl;
  std::cout << "  > Iterations: " << solver.iters() << std::endl;
  result -= ref_result;
  std::cout << "  > Relative deviation from result: " << viennacl::linalg::norm_2(result) / viennacl::linalg::norm_2(ref_result) << std::endl;
}


void run_amg(viennacl::linalg::bicgstab_tag & cg_solver,
             viennacl::vector<double> & vcl_vec,
             viennacl::vector<double> & vcl_result,
             viennacl::compressed_matrix<double> & vcl_compressed_matrix,
             viennacl::linalg::amg_tag & amg_tag)
{
  unsigned int coarselevels = amg_tag.get_coarselevels();
  

  viennacl::linalg::amg_precond<viennacl::compressed_matrix<double> > vcl_amg = viennacl::linalg::amg_precond<viennacl::compressed_matrix<double> > (vcl_compressed_matrix, amg_tag);
  std::cout << " * Setup phase (ViennaCL types)..." << std::endl;      
  vcl_amg.tag().set_coarselevels(coarselevels);
  vcl_amg.setup();
    
  std::cout << " * BiCGSTAB solver (ViennaCL types)..." << std::endl;         
  run_solver(vcl_compressed_matrix, vcl_vec, vcl_result, cg_solver, vcl_amg);

}







#if 0  
int main()
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
#endif
