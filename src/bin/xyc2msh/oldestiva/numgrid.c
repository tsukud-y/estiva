#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "estiva.h"
#include "confary.h"
#include "numgrid.h"
#include "op.h"

static int count=0;

void numgrid_initgamma(void)
{ count = 0;}

static int neq(char *str1, char *str2)
{
  if(str1 == (char *)NULL && str2 == (char *)NULL) return 0;
  if(str1 == (char *)NULL || str2 == (char *)NULL) return 1;
  return strcmp(str1,str2);
}
int numgrid_gamma(int *p,char *str)
{ 
  if(dim1(Mid) < ++count) return 0;

  while(neq(Mid[count].l,str))if(count<dim1(Mid)) count++;else return 0;
  *p = count;
  return 1;
}

void lambda(float **A, int i)
{ A[i][i]=1.0e+30;}

struct element  *Ele;
struct vertex   *Ver;
struct midpoint *Mid;

void readfile(char *name, char *name1, char *name2)
{ FILE *fp, *fp1, *fp2;

  if(cmpop("-h","40"))
    ary1(Ele,3200+1);  
  else if(cmpop("-h","20"))
    ary1(Ele,800+1);  
  else 
    ary1(Ele,200+1);  

  fprintf(stderr,"Ele = %d ",dim1(Ele));
    
  if(NULL == (fp1 = fopen(name1,"r"))){ perror(name1); exit(1);}
  { int i,e=dim1(Ele),a,b,c,ma,mb,mc;
    while(fscanf(fp1,"%d %d %d %d %d %d %d",&i,&a,&b,&c,&ma,&mb,&mc) != EOF)
      if(i<=e){
	Ele[i].a=a;Ele[i].b=b;Ele[i].c=c;
	Ele[i].ma=ma;Ele[i].mb=mb;Ele[i].mc=mc;
      }
  }fclose(fp1);

  if(NULL == (fp2 = fopen(name2,"r"))){ perror(name2); exit(1);}
  { int i,e=dim1(Ele),a,b,c,ua,ub,uc;
    while(fscanf(fp2,"%d %d %d %d %d %d %d",&i,&a,&b,&c,&ua,&ub,&uc) != EOF)
      if(i<=e){	
	Ele[i].ua=ua;Ele[i].ub=ub;Ele[i].uc=uc;
      }
  }fclose(fp2);
  { int e,v=0,m=0;
    forall(1,e,dim1(Ele)){
      v=more(v,Ele[e].a);v=more(v,Ele[e].b);v=more(v,Ele[e].c);
      m=more(m,Ele[e].ma);m=more(m,Ele[e].mb);m=more(m,Ele[e].mc);
    }
    ary1(Ver,v+1); ary1(Mid,m+1); 
    fprintf(stderr,"Ver = %d Mid = %d\n",dim1(Ver),dim1(Mid));
  }
  if(NULL == (fp = fopen(name,"r"))){ perror(name); exit(1);}
  { int i,v=dim1(Ver); float x, y;
    while(fscanf(fp,"%d%f%f",&i,&x,&y) != EOF)if(i<=v){Ver[i].x=x;Ver[i].y=y;}
  }fclose(fp);

  { static float xmin= -1.0, xmax=1.0, ymin= -1.0, ymax=1.0;
    static char G1[]="G1", G2[]="G2", G3[]="G3", G4[]="G4";
    int i;double x,y;
    forall(1,i,dim1(Ver)){
      x=Ver[i].x;y=Ver[i].y;
      Ver[i].l = NULL;
      if(x == xmin) Ver[i].l = G1;
      if(x == xmax) Ver[i].l = G2;
      if(y == ymin) Ver[i].l = G3;
      if(y == ymax) Ver[i].l = G4;
      if(x == xmin && y == ymin) Ver[i].l = NULL;
      if(x == xmin && y == ymax) Ver[i].l = NULL;
      if(x == xmax && y == ymin) Ver[i].l = NULL;
      if(x == xmax && y == ymax) Ver[i].l = NULL;
    }
  }
  { int e;double ax,ay,bx,by,cx,cy;
    forall(1,e,dim1(Ele)){
      ax=Ver[Ele[e].a].x;ay=Ver[Ele[e].a].y;
      bx=Ver[Ele[e].b].x;by=Ver[Ele[e].b].y;
      cx=Ver[Ele[e].c].x;cy=Ver[Ele[e].c].y;
      Ele[e].x= (ax+bx+cx)/3.0;
      Ele[e].y= (ay+by+cy)/3.0;
      Ele[e].s= ((ax*by+ay*cx+bx*cy)-(by*cx+ay*bx+cy*ax))/2.0;
    }
  }
  { int e,a,b,c,ua,ub,uc,ma,mb,mc;
    forall(1,e,dim1(Ele)){
      a=Ele[e].a;b=Ele[e].b;c=Ele[e].c;
      ua=Ele[e].ua;ub=Ele[e].ub;uc=Ele[e].uc;
      ma=Ele[e].ma;mb=Ele[e].mb;mc=Ele[e].mc;
      Mid[ma].l=NULL;Mid[mb].l=NULL;Mid[mc].l=NULL;
      if(ua==0){
	if(Ver[b].l != NULL) Mid[ma].l=Ver[b].l;
	if(Ver[c].l != NULL) Mid[ma].l=Ver[c].l;
      }
      if(ub==0){
	if(Ver[c].l != NULL) Mid[mb].l=Ver[c].l;
	if(Ver[a].l != NULL) Mid[mb].l=Ver[a].l;
      }
      if(uc==0){
	if(Ver[a].l != NULL) Mid[mc].l=Ver[a].l;
	if(Ver[b].l != NULL) Mid[mc].l=Ver[b].l;
      }
    }
  }
}

int numgrid_e;
#define val(n) forall(1,i,n) v[i]=va_arg(ap,double)
void FEM(void *v0,...)
{ double v[10]; va_list ap;
  int i,e;
  e=numgrid_e;
  if(e==1) reset(v0);
  va_start(ap,v0);

  if(dim2(v0)==dim1(Mid)&&dim1(v0)==dim1(Ele)){ 
    float **H;
    cp(v0,H); val(3);
    H[Ele[e].ma][e] += v[1];
    H[Ele[e].mb][e] += v[2];
    H[Ele[e].mc][e] += v[3];
  }
  else if(dim2(v0)==dim1(Mid)&&dim1(v0)==dim1(Mid)){
    float **K; int ma,mb,mc;
    cp(v0,K); val(9);
    ma=Ele[e].ma, mb=Ele[e].mb, mc=Ele[e].mc;

    K[ma][ma] += v[1], K[ma][mb] += v[2], K[ma][mc] += v[3];
    K[mb][ma] += v[4], K[mb][mb] += v[5], K[mb][mc] += v[6];
    K[mc][ma] += v[7], K[mc][mb] += v[8], K[mc][mc] += v[9];
  }
  else if(dim2(v0)== -1&&dim1(v0)==dim1(Ele)){
    float *r;
    cp(v0,r); val(1);
    r[e] += v[1];
  }
  else if(dim2(v0)== -1&&dim1(v0)==dim1(Mid)){
    float *M;
    cp(v0,M); val(9);
    M[Ele[e].ma] += v[1];
    M[Ele[e].mb] += v[5];
    M[Ele[e].mc] += v[9];
  }
  else if(dim2(v0)==dim1(Ele)&&dim1(v0)==dim1(Ele)){
    float **N; int ua,ub,uc;
    cp(v0,N); val(3);
    ua=Ele[e].ua;
    ub=Ele[e].ub;
    uc=Ele[e].uc;

    if(ua!=0) N[e][ua] += v[1];
    if(ub!=0) N[e][ub] += v[2];
    if(uc!=0) N[e][uc] += v[3];
  }

  va_end(ap);
}
#undef val
