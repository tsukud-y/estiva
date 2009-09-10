#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva.h"
#include "stack.h"
#include "foreach.h"
#include "confary.h"
#include "Delaunay.h"
#include "FILE.h"


#if 0
typedef struct { double x, y; char *label;}xyc;
typedef struct { int a,b,c,A,B,C,ma,mb,mc;}nde;
#endif

static xyc *Z;
static nde *N;


#define distance2(x0,y0,x1,y1) ((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))

#define sarrus(a11,a12,a13,a21,a22,a23,a31,a32,a33) \
 ((a11)*(a22)*(a33)+(a21)*(a32)*(a13)+(a31)*(a12)*(a23)\
 -(a13)*(a22)*(a31)-(a23)*(a32)*(a11)-(a33)*(a12)*(a21))

static double fplane(double x, double y, double z,
		     double x0,double y0,double z0,
		     double x1,double y1,double z1,
		     double x2,double y2,double z2)
{ double xn,yn,zn;

  xn = (y1-y0)*(z2-z0)-(z1-z0)*(y2-y0);
  yn = (z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);
  zn = (x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);

  return xn*x+yn*y+zn*z-(xn*x0+yn*y0+zn*z0);
}
static void cramer3(double *px,double *py,double *pz, 
		    double a11,double a12,double a13, 
		    double a21,double a22,double a23, 
		    double a31,double a32,double a33,
		    double b1, double b2, double b3 )
{ double det;

  *px = sarrus(b1 ,a12,a13, b2 ,a22,a23, b3 ,a32,a33);
  *py = sarrus(a11,b1 ,a13, a21,b2 ,a23, a31,b3 ,a33);
  *pz = sarrus(a11,a12,b1 , a21,a22,b2 , a31,a32,b3 );
  det = sarrus(a11,a12,a13, a21,a22,a23, a31,a32,a33);
  if(det != 0.0){ *px/=det;*py/=det;*pz/=det;}
}
static void rotate_right(int i)
{ int a,b,c,A,B,C;
  a=N[i].a; b=N[i].b; c=N[i].c; A=N[i].A; B=N[i].B; C=N[i].C;
  N[i].a=c; N[i].b=a; N[i].c=b; N[i].A=C; N[i].B=A; N[i].C=B;
}
static void rotate_left(int i)
{ int a,b,c,A,B,C;
  a=N[i].a; b=N[i].b; c=N[i].c; A=N[i].A; B=N[i].B; C=N[i].C;
  N[i].a=b; N[i].b=c; N[i].c=a; N[i].A=B; N[i].B=C; N[i].C=A;
}
static int incircle(int p,int e2)
{ double a, b, c, x,y,x0,y0,x1,y1,x2,y2;

  x=Z[p].x; y=Z[p].y;
  x0=Z[N[e2].a].x;  y0=Z[N[e2].a].y;
  x1=Z[N[e2].b].x;  y1=Z[N[e2].b].y;
  x2=Z[N[e2].c].x;  y2=Z[N[e2].c].y;

  cramer3(&a,&b,&c, 
	  x0,y0,1.0,
	  x1,y1,1.0,
	  x2,y2,1.0, 
	  x0*x0 + y0*y0,
	  x1*x1 + y1*y1, 
	  x2*x2 + y2*y2);

  return distance2(x0,y0,a/2.0,b/2.0)-distance2(x,y,a/2.0,b/2.0)>0.0? 1:0;
}
#define N_set(n,i,j,k,I,J,K) \
 N[n].a=i;N[n].b=j;N[n].c=k;N[n].A=I;N[n].B=J;N[n].C=K;


static int notLawson(double x, double y)
{
  int e,n;
  n = dim1(N);

  for(e=1;e<=n;e++){
    int a,b,c; double x0,y0,x1,y1,x2,y2;
    a = N[e].a; b = N[e].b; c = N[e].c; 
    x0 = Z[a].x; x1 = Z[b].x; x2 = Z[c].x;
    y0 = Z[a].y; y1 = Z[b].y; y2 = Z[c].y;

    if     (fplane(x,y,0.0, x0,y0,1.0, x1,y1,0.0, x2,y2,0.0)>0.0) continue;
    else if(fplane(x,y,0.0, x0,y0,0.0, x1,y1,1.0, x2,y2,0.0)>0.0) continue;
    else if(fplane(x,y,0.0, x0,y0,0.0, x1,y1,0.0, x2,y2,1.0)>0.0) continue;
    else return e;
  }
}

static int Lawson(int e, double x, double y)
{ 
  int n, count = 0;
  int a,b,c; double x0,y0,x1,y1,x2,y2;

  n = dim1(N);
  while(1){
    a = N[e].a; b = N[e].b; c = N[e].c; 
    x0 = Z[a].x; x1 = Z[b].x; x2 = Z[c].x;
    y0 = Z[a].y; y1 = Z[b].y; y2 = Z[c].y;

    if(n<count++) return notLawson(x,y);

    if     (fplane(x,y,0.0, x0,y0,1.0, x1,y1,0.0, x2,y2,0.0)>0.0) e = N[e].A;
    else if(fplane(x,y,0.0, x0,y0,0.0, x1,y1,1.0, x2,y2,0.0)>0.0) e = N[e].B;
    else if(fplane(x,y,0.0, x0,y0,0.0, x1,y1,0.0, x2,y2,1.0)>0.0) e = N[e].C;
    else return e;
  }
}
static int degeneracy(int e1,int e2)
{ int a, d, W, Y; double x1,y1,x2,y2,x3,y3,x4,y4;

  if(N[e1].B == e2) rotate_left(e1); if(N[e1].C == e2) rotate_right(e1);
  if(N[e2].B == e1) rotate_left(e2); if(N[e2].C == e1) rotate_right(e2);

  x1= Z[N[e1].a].x;  y1= Z[N[e1].a].y;
  x2= Z[N[e2].a].x;  y2= Z[N[e2].a].y;
  x3= Z[N[e2].b].x;  y3= Z[N[e2].b].y;
  x4= Z[N[e2].c].x;  y4= Z[N[e2].c].y;
 
  if(distance2(x1,y1,x2,y2)==distance2(x3,y3,x4,y4)) return 1;
  
  a=N[e1].a; Y=N[e1].C;
  d=N[e2].a; W=N[e2].C;
  N[e1].b=d; N[e1].A=W; N[e1].C=e2;
  N[e2].b=a; N[e2].A=Y; N[e2].C=e1;

  if(N[W].B == e2) rotate_left(W); if(N[W].C == e2) rotate_right(W);
  if(N[Y].B == e1) rotate_left(Y); if(N[Y].C == e1) rotate_right(Y);
  N[W].A = e1; N[Y].A = e2;

  return 0;
}
#define finv(g,f,n) ary1(g,n+1); forall(1,i,dim1(f)) g[f[i]]=i;

void estiva_Delaunay(xyc **Zo, nde **No)
     /* Delaunay(Z,N) estiva_Delaunay(&(Z),&(N)) */
{ int i,p,n,e,m,a,b,c,A,B,C,z; double xmin,ymin,xmax,ymax,length;
  static char super_node[] = "super_node";
  static int *fN,*fNinv,*fZ,*fZinv;
  Z= *Zo; n=1; z=dim1(Z); z -=3;

  
  xmin = Z[1].x; ymin = Z[1].y;
  xmax = Z[1].x; ymax = Z[1].y;
  forall(2,i,z){
    xmin=less(xmin,Z[i].x); xmax=more(xmax,Z[i].x);
    ymin=less(ymin,Z[i].y); ymax=more(ymax,Z[i].y);
  }



  length = more(xmax-xmin, ymax-ymin);
  Z[z+1].x = xmin - length/2.0;
  Z[z+1].y = ymin - length/2.0;
  Z[z+2].x = xmin + length*3.0;
  Z[z+2].y = ymin - length/2.0;
  Z[z+3].x = xmin - length/2.0;
  Z[z+3].y = ymin + length*3.0;
  Z[z+1].label = super_node;
  Z[z+2].label = super_node;
  Z[z+3].label = super_node;
  
  ary1(N,dim1(Z)*2);

  N_set(1,z+1,z+2,z+3,0,0,0);
  n=1;

  forall(1,i,z){
    int e0,e1,e2,a,b,c,A,B,C; 
    e0=Lawson(n,Z[i].x,Z[i].y); e1=n+1; e2=n+2; n+=2;
    a=N[e0].a,b=N[e0].b,c=N[e0].c,A=N[e0].A,B=N[e0].B,C=N[e0].C;
    
    N_set( 0,0,0,0, 0, 0, 0);
    N_set(e0,i,b,c,A ,e1,e2);  
    N_set(e1,a,i,c,e0,B, e2);
    N_set(e2,a,b,i,e0,e1, C);
    
    foreach(a)&N[B].A,&N[B].B,&N[B].C,end if(a==e0) a=e1;
    foreach(a)&N[C].A,&N[C].B,&N[C].C,end if(a==e0) a=e2;
    
    push(e0); push(e1); push(e2); 
    while(pop(e1))
      if(incircle(i,(e2=(N[e1].a==i?N[e1].A:(N[e1].b==i?N[e1].B:N[e1].C)))))
	if(!degeneracy(e1,e2)){ push(e1); push(e2);}
  }


  forall(1,i,dim1(N)){
    if(Z[N[i].a].label==NULL||Z[N[i].b].label==NULL||Z[N[i].c].label==NULL)
      continue;
    foreach(A)&N[i].A,&N[i].B,&N[i].C,end 
      foreach(a)&N[A].A,&N[A].B,&N[A].C,end if(a==i)a=0;
    foreach(a)&N[i].a,&N[i].b,&N[i].c,&N[i].A,&N[i].B,&N[i].C,end a=0;
  }

  n=0;forall(1,i,dim1(N))if(N[i].a!=0) n++;
  ary1(fN,n+1); n=0;forall(1,i,dim1(fN)){ while(N[n].a==0)n++; fN[i]=n++;}

  { static nde *N1;
    ary1(N1,dim1(fN)+1); forall(1,i,dim1(N1)) cp(N[fN[i]],N1[i]);
    ary1(N, dim1(fN)+1); forall(1,i,dim1(N))  cp(N1[i],N[i]);
  }

  /*
  finv(fNinv,fN,dim1(N)); 
  finv(g,f,n) 
  */

  ary1(fNinv,n+1); 
  if(fNinv == NULL) fprintf(stderr,"Can't malloc()\n");
  forall(1,i,dim1(fN)) fNinv[fN[i]]=i;
  forall(1,i,dim1(N)) foreach(A)&N[i].A,&N[i].B,&N[i].C,end A=fNinv[A];


  ary1(fN,0); 

  ary1(fNinv,0); /* error something wrong! */


  z=0; ary1(fZ,dim1(Z)+1); forall(1,i,dim1(N)) 
    foreach(a)&N[i].a,&N[i].b,&N[i].c,end if(fZ[a]==0)fZ[a]= ++z;
  forall(1,i,dim1(N)) foreach(a)&N[i].a,&N[i].b,&N[i].c,end a=fZ[a];
  finv(fZinv,fZ,z);
  { static xyc *Z1;
    ary1(Z1,z+1); forall(1,i,dim1(Z1)) cp(Z[fZinv[i]],Z1[i]);
    ary1(Z,z+1);  forall(1,i,dim1(Z))  cp(Z1[i],Z[i]);
  }
  ary1(fZ,0); ary1(fZinv,0);

  cp(Z,*Zo); cp(N,*No);
}
xyc *estiva_fp2Z(FILE *fp)
     /* fp2Z(fp) estiva_fp2Z(fp) */
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
#ifdef dela
#include "op.h"
static char *label(xyc *Z,int i)
{ static char *normal = "";
  return Z[i].label==NULL? normal:Z[i].label;
}
int main(int argc, char **argv)
{ FILE *fp; xyc *Z; nde *N; int i,z;
  initop(argc,argv);
  if(NULL==(fp=fopen(argv[1],"r"))) exit(1);
  
  Z = fp2Z(fp); fclose(fp);
  Delaunay(Z,N);
  
  if(defop("-o")){
    fp=stdout;
    fprintf(fp,"<xyc>\n");forall(1,i,dim1(Z))
      fprintf(fp,"%4d %f %f %s\n",i,Z[i].x,Z[i].y,label(Z,i));
    fprintf(fp,"\n<nde>\n");forall(1,i,dim1(N))
      fprintf(fp,"%4d  %4d %4d %4d  %4d %4d %4d\n",
	      i, N[i].a,N[i].b,N[i].c,N[i].A,N[i].B,N[i].C);
    fprintf(fp,"\n");
  }
  else{
    fp= popen("xgraph =+0+0","w");
    forall(1,i,dim1(N)) fprintf(fp,"move %f %f\n%f %f\n%f %f\n%f %f\n",
				Z[N[i].c].x,Z[N[i].c].y,	    
				Z[N[i].a].x,Z[N[i].a].y,
				Z[N[i].b].x,Z[N[i].b].y,
				Z[N[i].c].x,Z[N[i].c].y);
    pclose(fp);
  }
}
#endif
