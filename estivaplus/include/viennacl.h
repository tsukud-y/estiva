#ifndef _VIENNACL_H_
#define _VIENNACL_H_
#include <iostream>
#include <viennacl/scalar.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/gmres.hpp>
#include <viennacl/linalg/bicgstab.hpp>
#include <viennacl/linalg/compressed_matrix_operations.hpp>
#include <viennacl/linalg/row_scaling.hpp>
#include <viennacl/linalg/jacobi_precond.hpp>
#include <viennacl/linalg/amg.hpp>


using viennacl::linalg::amg_precond;
using viennacl::linalg::amg_tag;
using viennacl::linalg::ilut_precond;
using viennacl::linalg::jacobi_precond;
using viennacl::linalg::row_scaling;
using viennacl::compressed_matrix;
using viennacl::linalg::cg_tag;
using viennacl::linalg::gmres_tag;
using viennacl::linalg::bicgstab_tag;
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
typedef amg_precond< SparseMatrix > AMG;
typedef row_scaling< SparseMatrix > Scaling;
typedef jacobi_precond< SparseMatrix > Jacobi;

#endif
