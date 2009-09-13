#include <stdio.h>
#include <estiva/tmpfile.h>


static int tmpcounter = 0;
static char filename[1000][100];
static FILE *tmpfp[1000];

FILE *estiva_tmpopen(void)
{
  sprintf(filename[tmpcounter],"%d",tmpcounter);
  tmpfp[tmpcounter] = fopen(filename[tmpcounter],"w");
  tmpcounter++;
  return tmpfp[tmpcounter -1];
}

char *estiva_tmpname(FILE *fp)
{
  int i;
  for ( i = 1000; i>=0; i--) if ( fp == tmpfp[i]) break;
  return filename[i];
}

void estiva_tmpclose(FILE *fp)
{
  char cmd[100];
  sprintf(cmd,"rm %s",tmpname(fp));
  system(cmd);
  fclose(fp);
}
