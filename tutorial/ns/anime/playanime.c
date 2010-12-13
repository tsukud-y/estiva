#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main(int argc, char **argv) {
  FILE *pp, *fp;
  int c;  
  static long sequence=0;
  static char filename[400];

  if (argc == 2 ) sequence = atoi(argv[1]);

  pp = popen("gnuplot","w");
  
  while (1) {
    printf("sequencde = %ld\n",sequence);
    sprintf(filename,"%ld.gnuplot",sequence++);
    fp = fopen(filename,"r");
    if ( NULL == fp ) return 0;

    while ( EOF != (c = fgetc(fp)) ){
      fprintf(pp,"%c",c);
    }
    fclose(fp);
    fflush(pp);
    usleep(100);
  }
  return 0;
}
