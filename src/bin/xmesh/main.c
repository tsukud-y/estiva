#include <stdio.h>
#include <unistd.h>

void sleep_forever(void) {
  while(1) sleep(24*60*60);
}

FILE *tmpopen(char *name) {
  char filename[1024];
  sprintf(filename,"/tmp/%s",name);
  return fopen(filename,"w");
}

char *tmpname(char *name) {
  return name;
}


int main(int argc, char **argv) {
  FILE *pp, *xyc2msh;
  int c;

  if (argc == 2 ) {
    char command[1024];
    sprintf(command,"cat %s|xyc2msh|xmsh>/tmp/tmpfile.plt",argv[1]);
    xyc2msh = popen(command,"w");
  }
  else {
    xyc2msh = popen("xyc2msh|xmsh>/tmp/tmpfile.plt","w");
    while (EOF != (c = getchar())) putc(c,xyc2msh);
  }

  pclose(xyc2msh);



  pp = popen("gnuplot","w");
  fprintf(pp,"set nokey\n");
  fprintf(pp,"plot '/tmp/tmpfile.plt' w l\n");
  fflush(pp);

  sleep_forever();
  return 0;
}
