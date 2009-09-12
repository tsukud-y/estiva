#include <stdio.h>
#include <estiva/ary.h>
#include <estiva/foreach.h>
#include <estiva/mesh.h>



static void makeMid(nde *N)
{
  int i, e, m, n, u;
  e = dim1(N);

  for(i=1;i<=e;i++)
    foreach(m)&N[i].A,&N[i].B,&N[i].C,end
      m = -m;

  n=1;
  for(i=1;i<=e;i++) 
    foreach(m)&N[i].A,&N[i].B,&N[i].C,end 
      if(m<=0){
	foreach(u)&N[-m].A,&N[-m].B,&N[-m].C,end if(u== -i) u=n;
	m=n++;
      }
}

#define max(x,y) (x>y?x:y)

static char *bound(char *s1, char *s2)
{
  static char *blank = "";

  if(s1 == NULL || s2 == NULL) return blank;
  
  return strcmp(s1,s2)<0?s1:s2;
}

static xyc *makeMV(xyc *Z, nde *N)
{
  static xyc *MV;
  int i, e, m, k, v;
  e = dim1(N);
  
  m = 1;
  for(i=1;i<=e;i++)
    foreach(k)&N[i].A,&N[i].B,&N[i].C,end
      m = max(m,k);

  ary1(MV,m+1);


  for(i=1;i<=e;i++){
    MV[N[i].A].label = bound(Z[N[i].b].label,Z[N[i].c].label);
    MV[N[i].B].label = bound(Z[N[i].c].label,Z[N[i].a].label);
    MV[N[i].C].label = bound(Z[N[i].a].label,Z[N[i].b].label);

    MV[N[i].A].x = (Z[N[i].b].x + Z[N[i].c].x)/2.0;
    MV[N[i].B].x = (Z[N[i].c].x + Z[N[i].a].x)/2.0;
    MV[N[i].C].x = (Z[N[i].a].x + Z[N[i].b].x)/2.0;

    MV[N[i].A].y = (Z[N[i].b].y + Z[N[i].c].y)/2.0;
    MV[N[i].B].y = (Z[N[i].c].y + Z[N[i].a].y)/2.0;
    MV[N[i].C].y = (Z[N[i].a].y + Z[N[i].b].y)/2.0;
  }

  return MV;
}

static xyc *makeM(xyc *Z,nde *N)
{
  makeMid(N);
  return makeMV(Z,N);
}

void* Ver2Mid(xyc* Z, nde * N)
{
  static xyc *Mid;
  Mid = makeM(Z,N);
  return Mid;
}
