#include <stdio.h>

static void cp_fp(FILE* in, FILE* out)
{
  int c;

  while(EOF !=(c = getc(in))) putc(c,out);
  fflush(out);
  rewind(out);
}


void *argf(int argc,char **argv)
{
  FILE *fp;

  for(; argc != 1; argc--)
    if(NULL != (fp=fopen(argv[argc-1],"r")) )
      return fp;

  fp = tmpfile();
  cp_fp(stdin, fp);

  return fp;
}
