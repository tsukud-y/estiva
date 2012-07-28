#include "Delaunay.h"
#undef push
#undef pop
#include <stack>
#include <cmath>
#include <unistd.h>



#define distance2(x0,y0,x1,y1) ((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))


static void rotate_right(int i,vector<Nde>&N)
{ int a,b,c,A,B,C;
  a=N[i].a; b=N[i].b; c=N[i].c; A=N[i].A; B=N[i].B; C=N[i].C;
  N[i].a=c; N[i].b=a; N[i].c=b; N[i].A=C; N[i].B=A; N[i].C=B;
}
static void rotate_left(int i,vector<Nde>&N)
{ int a,b,c,A,B,C;
  a=N[i].a; b=N[i].b; c=N[i].c; A=N[i].A; B=N[i].B; C=N[i].C;
  N[i].a=b; N[i].b=c; N[i].c=a; N[i].A=B; N[i].B=C; N[i].C=A;
}


static int degeneracy(int e1,int e2,vector<Xyc>&Z, vector<Nde>&N)
{ int a, d, W, Y; double x1,y1,x2,y2,x3,y3,x4,y4;

  if(N[e1].B == e2) rotate_left(e1,N); if(N[e1].C == e2) rotate_right(e1,N);
  if(N[e2].B == e1) rotate_left(e2,N); if(N[e2].C == e1) rotate_right(e2,N);

  x1= Z[N[e1].a].x;  y1= Z[N[e1].a].y;
  x2= Z[N[e2].a].x;  y2= Z[N[e2].a].y;
  x3= Z[N[e2].b].x;  y3= Z[N[e2].b].y;
  x4= Z[N[e2].c].x;  y4= Z[N[e2].c].y;
 
  if(distance2(x1,y1,x2,y2)==distance2(x3,y3,x4,y4)) return 1;
  
  a=N[e1].a; Y=N[e1].C;
  d=N[e2].a; W=N[e2].C;
  N[e1].b=d; N[e1].A=W; N[e1].C=e2;
  N[e2].b=a; N[e2].A=Y; N[e2].C=e1;

  if(N[W].B == e2) rotate_left(W,N); if(N[W].C == e2) rotate_right(W,N);
  if(N[Y].B == e1) rotate_left(Y,N); if(N[Y].C == e1) rotate_right(Y,N);
  N[W].A = e1; N[Y].A = e2;

  return 0;
}




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



static int incircle(int p,int e2,vector<Xyc>&Z,vector<Nde>&N)
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

void SplitT(long i, long e0, vector<Nde>&N)
{
  long e1, e2;
  e1 = N.size();
  e2 = N.size()+1;

  Nde nde;

  N.push_back(nde);
  N.push_back(nde);


  long a,b,c,A,B,C;


  a=N[e0].a; b=N[e0].b; c=N[e0].c; A=N[e0].A; B=N[e0].B; C=N[e0].C;   
  

#define N_set(n,i,j,k,I,J,K) \
  N[n].a=i;N[n].b=j;N[n].c=k;N[n].A=I;N[n].B=J;N[n].C=K;

  N_set( 0,0,0,0, 0, 0, 0);
  N_set(e0,i,b,c,A ,e1,e2);  
  N_set(e1,a,i,c,e0,B, e2);
  N_set(e2,a,b,i,e0,e1, C);

  if(N[B].A==e0) N[B].A=e1;
  if(N[B].B==e0) N[B].B=e1;
  if(N[B].C==e0) N[B].C=e1;

  if(N[C].A==e0) N[C].A=e2;
  if(N[C].B==e0) N[C].B=e2;
  if(N[C].C==e0) N[C].C=e2;
}


int main(int argc, char ** argv)
{
  vector<Xyc> Z;
  Xyc z;

  z.x=0.0, z.y=0.0, z.label="super_node"; Z.push_back(z);

#if 0
  z.x=0.0, z.y=0.0, z.label="G1"; Z.push_back(z);
  z.x=0.0, z.y=0.5, z.label="G1"; Z.push_back(z);
  z.x=0.5, z.y=0.0, z.label="G1"; Z.push_back(z);
  z.x=1.0, z.y=0.0, z.label="G1"; Z.push_back(z);
  z.x=1.0, z.y=0.5, z.label="G1"; Z.push_back(z);
  z.x=0.5, z.y=0.5, z.label="";   Z.push_back(z);
  z.x=0.0, z.y=1.0, z.label="G1"; Z.push_back(z);
  z.x=0.5, z.y=1.0, z.label="G1"; Z.push_back(z);
  z.x=1.0, z.y=1.0, z.label="G1"; Z.push_back(z);
#endif

  
  for ( z.x = 0.0; z.x < 1.0; z.x+=0.1)
    for ( z.y = 0.0; z.y < 1.0; z.y+=0.1)
      {
	z.label="";
	Z.push_back(z);
      }

  z.x=0.0, z.y=0.0, z.label="super_node"; Z.push_back(z);
  z.x=0.0, z.y=0.0, z.label="super_node"; Z.push_back(z);
  z.x=0.0, z.y=0.0, z.label="super_node"; Z.push_back(z);


  unsigned long n = Z.size();

  Z[n-3].x=Xmin(Z), Z[n-3].y=Ymax(Z);  Z[n-2].x=Xmax(Z), Z[n-2].y=Ymax(Z);
  Z[  0].x=Xmin(Z), Z[  0].y=Ymin(Z);  Z[n-1].x=Xmax(Z), Z[n-1].y=Ymin(Z);

  double length = sqrt( pow(Xmax(Z)-Xmin(Z),2) + pow(Ymax(Z)-Ymin(Z),2))/2;

  Z[n-3].x-=length, Z[n-3].y+=length;  Z[n-2].x+=length, Z[n-2].y+=length;
  Z[  0].x-=length, Z[  0].y-=length;  Z[n-1].x+=length, Z[n-1].y-=length;

  
  vector<Nde> N;
  Nde nde;
 
  nde.a=  0, nde.b=0,   nde.c=0,   nde.A=0, nde.B=0, nde.C=0;
  N.push_back(nde);
  nde.a=  0, nde.b=n-1, nde.c=n-3, nde.A=2, nde.B=0, nde.C=0;
  N.push_back(nde);
  nde.a=n-2, nde.b=n-3, nde.c=n-1, nde.A=1, nde.B=0, nde.C=0;
  N.push_back(nde);




  FILE *pp = popen("gnuplot","w");

  

  long i = 1;
  for ( i = 1 ; i <Z.size()-3; i++) {
  long e0, e1,e2,e3;

  e0 = SearchT(Z,N,i);
  SplitT(i,e0,N);

  Putmesh(Z,N);
  Xmesh(pp,Z,N);
  sleep(1);

  using namespace std;
  stack<long> st;

  
  if ( 0 != N[e0].A ) st.push(N[e0].A);
  if ( 0 != N[e0].B ) st.push(N[e0].B);
  if ( 0 != N[e0].C ) st.push(N[e0].C);
  st.push(e0);

  
  while(!st.empty()){
    e1 = st.top(); st.pop();
    while (e1 == 0 ) {e1=st.top(); st.pop();}

    if ( N[e1].a == i ) e2 = N[e1].A;
    if ( N[e1].b == i ) e2 = N[e1].B;
    if ( N[e1].c == i ) e2 = N[e1].C;

    if (e2 != 0 )
      if(incircle(i,e2,Z,N))
	if(!degeneracy(e1,e2,Z,N)){ st.push(e1); st.push(e2);}
    N_set( 0,0,0,0, 0, 0, 0);
  }

  
  Putmesh(Z,N);
  Xmesh(pp,Z,N);
  sleep(1);
  }

  Putmesh(Z,N);
  Xmesh(pp,Z,N);
  sleep(300);
}
