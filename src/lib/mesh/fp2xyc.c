#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva/std.h"
#include <estiva/ary.h>
#include "estiva/mesh.h"

#define  forFILE(fp) while(estiva_forFILE(fp)) 
#define  S(n) estiva_S(n) 
#define  FILE_cp(fp,out) estiva_FILE_cp(fp,out) 

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
static int estiva_forFILE(FILE *fp)
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
static char *estiva_S(int n)
     /* S(n) estiva_S(n) */
{ 
  if(n==1) return S1[0]=='\0'?NULL:S1;
  if(n==2) return S2[0]=='\0'?NULL:S2;
  if(n==3) return S3[0]=='\0'?NULL:S3;
  if(n==0) return S0;
  return NULL;
}
static void estiva_FILE_cp(FILE *fp,FILE *out)
     /* FILE_cp(fp,out) estiva_FILE_cp(fp,out) */
{ char c; 
  while(EOF !=(c=getc(fp))) putc(c,out);
}     

xyc *estiva_fp2xyc(FILE *fp)
/* fp2xyc(fp) estiva_fp2xyc(fp) */
{ static xyc *Z; int i,z; FILE *tfp;
  
  FILE_cp(fp,(tfp=tmpfile())); rewind(tfp);
  z=0; forFILE(tfp)if(S(1)!=NULL&&S(2)!=NULL)z++; rewind(tfp);
  
  ary1(Z,z+4);
  i= 1; forFILE(tfp)if(S(1)!=NULL&&S(2)!=NULL){
    Z[i].x= atof(S(1)); Z[i].y= atof(S(2));
    Z[i].label= (S(3)==NULL?NULL:strdup(S(3)));
    i++;
  }fclose(tfp);
  return Z;
}
