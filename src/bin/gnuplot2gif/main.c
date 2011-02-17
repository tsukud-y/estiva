#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

char *num2str(long n){
  long i;
  static char ret[100];
  sprintf(ret,"00000000%ld",n);
  i = strlen( ret );
  
  return &ret[i-8];
}


int main(int argc, char **argv) {
  FILE *pp, *fp;
  int c;  
  static long sequence=0;
  static char filename[400], cmd[400];

  if (argc == 2 ) sequence = atoi(argv[1]);

  pp = popen("gnuplot","w");
  
  while (1) {
    printf("sequencde = %ld\n",sequence);
    sprintf(filename,"%s.gnuplot",num2str(sequence++));
    fp = fopen(filename,"r");
    if ( NULL == fp ) return 0;

    fprintf(pp,"set terminal gif\n");
    sprintf(cmd,"set output \"%s.gif\"\n",num2str(sequence));
    fprintf(pp,"%s",cmd);
    while ( EOF != (c = fgetc(fp)) ){
      fprintf(pp,"%c",c);
    }
    fclose(fp);
    fflush(pp);
    usleep(100);
  }
  return 0;
}
