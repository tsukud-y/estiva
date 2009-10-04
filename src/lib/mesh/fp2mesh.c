#include <stdio.h>
#include <stdlib.h>
#include <estiva/ary.h>
#include <estiva/mesh.h>

static void cp_fp(FILE* in, FILE* out)
{
  int c;
  
  while(EOF !=(c = getc(in))) putc(c,out);
  fflush(out);
  rewind(out);
}

void fp2mesh(FILE* fp, xyc** Zp, nde** Np)
{
  static xyc* Z; static nde* N;
  char buf1000[1000], buf200[200];
  int i, j, m=0, n=0, state = 1;
  

  while(fgets(buf1000,999,fp)) switch(state){

  case 1:
    sscanf(buf1000,"%s",buf200);
    if(!strcmp(buf200,"<nde>")) state = 2;
    else  m = atoi(buf200);
    break;

  case 2:
    sscanf(buf1000,"%s",buf200);
    n = atoi(buf200);
    break;
  }

  ary1(Z,m+1); ary1(N,n+1);
  
  if(Z == NULL || N == NULL){
    fprintf(stderr,"poisson: Can't alloc memory!\n");
    abort();
  }

  rewind(fp);
  
  state = 1;

  while(fgets(buf1000,999,fp)) switch(state){
    
  case 1:
    sscanf(buf1000,"%s",buf200);
    if(!strcmp(buf200,"<xyc>")) state = 1;
    else if(!strcmp(buf200,"<nde>")) state = 2;
    else{
      sscanf(buf200,"%d",&i);
      Z[i].label = (char*)malloc(16);
      if(Z[i].label == NULL){
	fprintf(stderr,"poisson: Can't alloc memory\n");
	abort();
      }
      sscanf(buf1000,"%d %lf %lf %s",&m, &Z[i].x, &Z[i].y,Z[i].label);
    }
    break;

  case 2:
    sscanf(buf1000,"%s",buf200);
    sscanf(buf200,"%d",&j);
    sscanf(buf1000,"%d %d %d %d %d %d %d",&n, &N[j].a, &N[j].b, &N[j].c,
	   &N[j].A,&N[j].B,&N[j].C);
    break;
  }
  *Zp = Z;
  *Np = N;
}
