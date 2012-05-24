using namespace std;

extern "C" {
  extern void stwart_(void*,void*,void*,void*,void*,void*,void*,void*,
                      void*,void*,void*,void*,void*,void*,void*,void*,void*);
};
extern void Solverorg(Matrix&,Vector&,Vector&);
extern Matrix makeA(void);
extern void printMatrix(Matrix&);
extern void printVector(Vector&);
extern void use(long);
extern Vector makeb(long);
extern Matrix TransMatrix(Matrix &);
extern std::vector<long> genMA(Matrix &);
extern std::vector<long> genIA(Matrix &, std::vector<long>&);
extern void Stwart(std::vector<long>&,std::vector<long>&,
		   std::vector<long>&,std::vector<long>&);
extern Vector genb0(Vector &, std::vector<long> &);
extern Matrix genA0(Matrix &, std::vector<long>&, std::vector<long>&);
extern Vector genx(Vector &, std::vector<long>&);
extern void StwartMethod(Matrix &,std::vector<long>&,std::vector<long>&);
extern void SolverWithStwart(Matrix &, Vector &, Vector &);

static map<unsigned int, double>::iterator itr;

#define forMatrix(A,i,j)						\
  for (unsigned long j, i=0; i<A.size(); i++)				\
    for (itr=A[i].begin(), j=0; j=itr->first, itr != A[i].end(); itr++) \
      if ( A[i][j] != 0.0 )

#define forVector(b,i)  for ( unsigned long i = 0; i < b.size(); i++)
