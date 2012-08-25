#include <iostream>
#include <estivaplus.h>

using namespace std;

struct Xyc { double x; double y; string label;};
struct Nde { long a,b,c,A,B,C;};

class Mesh {
 private:
  static void SplitT(long i, long e0, vector<Nde>&N);
  static long SearchT(vector<Xyc>&Z, vector<Nde>&N, long i);
  static int incircle(int p,int e2,vector<Xyc>&Z,vector<Nde>&N);
  static int degeneracy(int e1,int e2,vector<Xyc>&Z, vector<Nde>&N);
  static void GenSuperNodes(vector<Xyc>&Z, vector<Nde>&N);
  static void DelaunayAlgo(vector<Xyc>&Z,vector<Nde>&N);
  static void VanishSuperNodes(vector<Xyc>&Z,vector<Nde>&N);
  static void VanishBT(vector<Xyc>&Z,vector<Nde>&N);
  static void SortTri(vector<Xyc> Z, vector<Nde>& N);
  static void GenRelation(vector<Xyc> Z, vector<Nde>&N);
  static void Normalization(vector<Xyc>&Z,vector<Nde>&N);
  static void Polynomial2(vector<Xyc>&Z, vector<Nde>&N);
  static double Xmin(vector<Xyc>&Z);
  static double Xmax(vector<Xyc>&Z);
  static double Ymin(vector<Xyc>&Z);
  static double Ymax(vector<Xyc>&Z);
  static int incircleDelta(int p,int e2,vector<Xyc>&Z,vector<Nde>&N);
  static void XMeshP2(FILE *pp,vector<Xyc>&Z, vector<Nde>&N);
 public:
  static void Gen(vector<Xyc>&Z, vector<Nde>&N);
  static void FPut(FILE *fp, vector<Xyc>&Z, vector<Nde>&N);
  static void X(FILE *fp, vector<Xyc>&Z, vector<Nde>&N);
};

