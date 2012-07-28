#include <estivaplus.h>

using namespace std;

struct Xyc { double x; double y; string label;};
struct Nde { long a,b,c,A,B,C;};

double Xmin(vector<Xyc>&Z);
double Xmax(vector<Xyc>&Z);
double Ymin(vector<Xyc>&Z);
double Ymax(vector<Xyc>&Z);

void Putmesh(vector<Xyc>&Z, vector<Nde>&N);
void Xmesh(FILE *fp, vector<Xyc>&Z, vector<Nde>&N);
long SearchT(vector<Xyc>&Z, vector<Nde>&N, long i);
