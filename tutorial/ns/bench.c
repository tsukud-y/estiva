#include <stdio.h>


int main(){
  char cmd[1000];
  long n;

  
  for ( n = 1; n <=200; n++){
    sprintf(cmd,
	    "time ./a.out -plotscale 0.001 -kn 1 -n %ld 2>foo",n); 
    system(cmd);
    system("cat foo >> bench.txt");
  }

}
