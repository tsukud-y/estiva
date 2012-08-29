#include "estivaplus/Mesh.h"

int main(int argc, char ** argv)
{
  initop(argc,argv);
  vector<Xyc> Z; vector<Nde> N;  

  Mesh::Rect(Z,N);
  //Mesh::FPut(stdout,Z,N);
  Mesh::Gnuplot(Z,N);
}
