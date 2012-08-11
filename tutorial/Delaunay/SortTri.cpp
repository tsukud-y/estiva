#include "Delaunay.h"
#undef push
#undef pop
#include <stack>
#include <cmath>
#include <unistd.h>
#include <algorithm>

using namespace std;

static vector<Xyc> Zv;

class LessInt {
public:
  bool operator()(const Nde& riLeft, const Nde& riRight) const {
    Xyc GLeft, GRight;
    
    GLeft.x = (Zv[riLeft.a].x + Zv[riLeft.b].x + Zv[riLeft.c].x) / 3;
    GLeft.y = (Zv[riLeft.a].y + Zv[riLeft.b].y + Zv[riLeft.c].y) / 3;

    GRight.x = (Zv[riRight.a].x + Zv[riRight.b].x + Zv[riRight.c].x) / 3;
    GRight.y = (Zv[riRight.a].y + Zv[riRight.b].y + Zv[riRight.c].y) / 3;
    
    if ( GLeft.y == GRight.y )
      return GLeft.x < GRight.x;
    return GLeft.y < GRight.y;
  }
};

void SortTri(vector<Xyc> Z, vector<Nde>& N)
{
  Zv = Z;
  sort(N.begin()+1,N.end(),LessInt());
}


