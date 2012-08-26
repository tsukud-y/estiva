#include "estivaplus/Mesh.h"

void Mesh::Gnuplot(vector<Xyc>&Z, vector<Nde>&N)
{
  FILE *pp = popen("gnuplot","w");
  if ( fork() == 0 ){
    Mesh::X(pp,Z,N);
    sleep(300);
  }
}
