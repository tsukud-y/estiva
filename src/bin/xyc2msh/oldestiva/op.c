#include <string.h>
#include "op.h"

static int argc;
static char **argv;

void op_initop(int pargc, char **pargv)
{
  argc = pargc; argv = pargv;
}

int op_defop(char *str)
{ static int i;
  for(i=argc-1; 0<i; i--)if(!strcmp(str,argv[i])) return 1;
  return 0;
}

char *op_getop(char *str)
{ static int i; 
  static char nullstr[] = "\0";

  for(i=argc-2;0<i;i--)if(!strcmp(str,argv[i])) return argv[i+1];

  if(!strcmp(str,"") && 1<argc) return argv[i+1];

  return nullstr;
}

int op_cmpop(char *opt, char *str)
{
  if(!strcmp(getop(opt),str)) return 1;
  return 0;
}

