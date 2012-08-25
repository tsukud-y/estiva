#include "Mesh.h"

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



int Mesh::degeneracy(int e1,int e2,vector<Xyc>&Z, vector<Nde>&N)
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
