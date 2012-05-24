#include "estivaplus.h"
#include "viennacl.h"
#include "stwart.h"

int main(){
  Matrix A  = makeA();
  Vector b  = makeb(A.size());
  Vector x(A.size());


  SolverWithStwart(A,x,b);

  b = viennacl::linalg::prod(A,x);
  printMatrix(A);
  printVector(b);
  return 0;
}
