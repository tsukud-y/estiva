#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva.h"
#include "confary.h"

#define reary1(a,n) ary1(a,n)
static char *S0;

static int getline(char *buf,FILE *fp)
{ char c; int i;
  for(i=0;1;){
    c = getc(fp);
    buf[i] = '\0';
    if(c=='\n'){ return 0;}
    if(c==EOF ){ return 1;}
    buf[i]=c;
    i++;
  }
}
static char S1[999], S2[999], S3[999];
int estiva_forFILE(FILE *fp)
     /* forFILE(fp) while(estiva_forFILE(fp)) */
{ int sup;
  sup=128; ary1(S0,sup);
  if(getline(S0,fp))if(feof(fp)){ary1(S0,0); return 0;}
  S1[0] = '\0';
  S2[0] = '\0';
  S3[0] = '\0';
  sscanf(S0,"%s %s %s",S1,S2,S3);
  return 1;
}
char *estiva_S(int n)
     /* S(n) estiva_S(n) */
{ 
  if(n==1) return S1[0]=='\0'?NULL:S1;
  if(n==2) return S2[0]=='\0'?NULL:S2;
  if(n==3) return S3[0]=='\0'?NULL:S3;
  if(n==0) return S0;
}
void estiva_FILE_cp(FILE *fp,FILE *out)
     /* FILE_cp(fp,out) estiva_FILE_cp(fp,out) */
{ char c; 
  while(EOF !=(c=getc(fp))) putc(c,out);
}     
