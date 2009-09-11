#include <stdio.h>
#include <estiva/ary.h>
#include <estiva/mesh.h>



main(){
  int i;
  static xyc *Z; static nde *N;
  
  fp2mesh(fopen("foo.mesh","r"),&Z,&N);


  for(i=1; i<=dim1(Z); i++)
    printf("%d %f %f %s\n",i,Z[i].x,Z[i].y, Z[i].label);

  for(i=1; i<=dim1(N); i++)
    printf("%d %d %d %d\n",i,N[i].a,N[i].b,N[i].c);
  
}
