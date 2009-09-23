#include <stdio.h>
#include <estiva/tmpfile.h>


static int tmpcounter = 0;
static char filename[100][100];
static FILE *tmpfp[100];

FILE *estiva_tmpopen(void)
{
  int i;
  for ( i = 0; i<=98; i++) if ( tmpfp[i] == NULL) break;
  sprintf(filename[i],"%d",tmpcounter);
  tmpfp[i] = fopen(filename[i],"w");
  tmpcounter++;
  return tmpfp[i];
}

char *estiva_tmpname(FILE *fp)
{
  int i;
  for ( i = 0; i<=98; i++) if ( fp == tmpfp[i]) break;
  return filename[i];
}

void estiva_tmpclose(FILE *fp)
{
  int i;
  char cmd[100];
  sprintf(cmd,"rm %s",tmpname(fp));
  system(cmd);
  fclose(fp);
  for ( i = 0; i<=98; i++) if ( fp == tmpfp[i]) tmpfp[i] = NULL;
}
