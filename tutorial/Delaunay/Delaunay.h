#include <iostream>
#include <estivaplus.h>

using namespace std;

struct Xyc { double x; double y; string label;};
struct Nde { long a,b,c,A,B,C;};

double Xmin(vector<Xyc>&Z);
double Xmax(vector<Xyc>&Z);
double Ymin(vector<Xyc>&Z);
double Ymax(vector<Xyc>&Z);

void FPutMesh(FILE *fp, vector<Xyc>&Z, vector<Nde>&N);
void XMesh(FILE *fp, vector<Xyc>&Z, vector<Nde>&N);
long SearchT(vector<Xyc>&Z, vector<Nde>&N, long i);
void GenSuperNodes(vector<Xyc>&Z, vector<Nde>&N);
void SplitT(long i, long e0, vector<Nde>&N);
int incircle(int p,int e2,vector<Xyc>&Z,vector<Nde>&N);
int incircleDelta(int p,int e2,vector<Xyc>&Z,vector<Nde>&N);
int degeneracy(int e1,int e2,vector<Xyc>&Z, vector<Nde>&N);
void VanishSuperNodes(vector<Xyc>&Z,vector<Nde>&N);
void DelaunayAlgo(vector<Xyc>&Z,vector<Nde>&N);
void VanishBT(vector<Xyc>&Z,vector<Nde>&N);
void SortTri(vector<Xyc> Z, vector<Nde>& N);
void GenRelation(vector<Xyc> Z, vector<Nde>&N);
void Normalization(vector<Xyc>&Z,vector<Nde>&N);
void Polynomial2(vector<Xyc>&Z, vector<Nde>&N);
void GenMesh(vector<Xyc>&Z, vector<Nde>&N);
void XMeshP2(FILE *pp,vector<Xyc>&Z, vector<Nde>&N);

#define distance2(x0,y0,x1,y1) ((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))
