#include <stdio.h>
#include <stdio.h>
#include <string.h>
#include <estiva/op.h>

static int argc;
static char **argv;

int    *estiva_pargc;
char ***estiva_pargv;

void estiva_initop(int *pargc, char ***pargv)
{
  argc = *pargc; argv = *pargv;
  estiva_pargc = pargc;
  estiva_pargv = pargv;
}

int estiva_defop(char *str)
{ static int i;
  for(i=argc-1; 0<i; i--)if(!strcmp(str,argv[i])) return 1;
  return 0;
}

char *estiva_getop(char *str)
{ static int i; 
  static char nullstr[] = "\0";

  for(i=argc-2;0<i;i--)if(!strcmp(str,argv[i])) return argv[i+1];

  if(!strcmp(str,"") && 1<argc) return argv[i+1];

  return nullstr;
}



static void cp_fp(FILE* in, FILE* out)
{
  int c;

  while(EOF !=(c = getc(in))) putc(c,out);
  fflush(out);
  rewind(out);
}


void *estiva_stdfp(void)
{
  FILE *fp;
  int i = argc;

  if (NULL != (fp=fopen(getop("-f"),"r")) )
    return fp;

  for(; i != 1; i--)
    if(NULL != (fp=fopen(argv[i-1],"r")) )
      return fp;

  fp = tmpfile();
  cp_fp(stdin, fp);

  return fp;
}
