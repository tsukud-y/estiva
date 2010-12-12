#include <stdio.h>


int main() {
  FILE *pp, *fp;
  int c;  
  static long sequence=0;
  static char filename[400];

  pp = popen("gnuplot","w");
  
  while (1) {
    sprintf(filename,"%ld.gnuplot",sequence++);
    fp = fopen(filename,"r");
    if ( NULL == fp ) return 0;

    while ( EOF != (c = fgetc(fp)) ){
      fprintf(pp,"%c",c);
    }
    fclose(fp);
    fflush(pp);
    sleep(1);
  }
  return 0;
}
