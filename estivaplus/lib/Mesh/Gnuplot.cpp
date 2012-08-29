#include "estivaplus/Mesh.h"

void Mesh::Gnuplot(vector<Xyc>&Z, vector<Nde>&N)
{
  if ( fork() == 0 ){
    FILE *pp = popen("gnuplot","w");
    Mesh::X(pp,Z,N);
    sleep(300);
    pclose(pp);
  }
}
