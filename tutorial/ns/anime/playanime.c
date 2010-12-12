#include <stdio.h>
#include <unistd.h>

int main() {
  FILE *pp, *fp;
  int c;  
  static long sequence=0;
  static char filename[400];

  pp = popen("gnuplot","w");
  
  while (1) {
    static long k = 0;
    sprintf(filename,"%ld.gnuplot",sequence++);
    fp = fopen(filename,"r");
    if ( NULL == fp ) return 0;

    while ( EOF != (c = fgetc(fp)) ){
      fprintf(pp,"%c",c);
    }
    fclose(fp);
    fflush(pp);
    printf("k = %ld\n",k++);
    fflush(stdout);
    usleep(100);
  }
  return 0;
}
