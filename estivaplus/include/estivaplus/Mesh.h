#include <iostream>
#include <cmath>
#include <algorithm>
#include <unistd.h>
#include <estivaplus.h>

using namespace std;

struct Xyc { double x; double y; string label;};
struct Nde { long a,b,c,A,B,C;};

class Mesh {
 private:
  static void SplitT(long i, long e0, vector<Nde>&N);
  static int incircle(int p,int e2,vector<Xyc>&Z,vector<Nde>&N);
  static int degeneracy(int e1,int e2,vector<Xyc>&Z, vector<Nde>&N);
  static long SearchT(vector<Xyc>&Z, vector<Nde>&N, long i);
  static void SortTri(vector<Xyc>&Z, vector<Nde>& N);
  static void VanishBT(vector<Xyc>&Z,vector<Nde>&N);
  static void GenRelation(vector<Xyc>&Z, vector<Nde>&N);
  static void Polynomial2(vector<Xyc>&Z, vector<Nde>&N);
  static void DelaunayAlgo(vector<Xyc>&Z,vector<Nde>&N);
  static void GenSuperNodes(vector<Xyc>&Z, vector<Nde>&N);
  static void Normalization(vector<Xyc>&Z, vector<Nde>&N);
  static void VanishSuperNodes(vector<Xyc>&Z, vector<Nde>&N);
  static void XP2(FILE *pp,vector<Xyc>&Z, vector<Nde>&N);
  static int incircleDelta(int p,int e2,vector<Xyc>&Z,vector<Nde>&N);
  static double Xmin(vector<Xyc>&Z);
  static double Xmax(vector<Xyc>&Z);
  static double Ymin(vector<Xyc>&Z);
  static double Ymax(vector<Xyc>&Z);
 public:
  static void Gen(vector<Xyc>&Z, vector<Nde>&N);
  static void FPut(FILE *fp, vector<Xyc>&Z, vector<Nde>&N);
  static void X(FILE *fp, vector<Xyc>&Z, vector<Nde>&N);
};

double distance2(double x0, double y0, double x1, double y1);
double sarrus(double a11,double a12,double a13,
	      double a21,double a22,double a23,
	      double a31,double a32,double a33);
void cramer3(double *px,double *py,double *pz,
	     double a11,double a12,double a13,
	     double a21,double a22,double a23,
	     double a31,double a32,double a33,
	     double b1, double b2, double b3 );
