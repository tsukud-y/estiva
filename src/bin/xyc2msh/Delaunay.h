#ifndef DELAUNAY_H
#define DELAUNAY_H
typedef struct { double x, y; char *label;}xyc;
typedef struct { int a,b,c,A,B,C,ma,mb,mc;}nde;
extern void estiva_Delaunay(xyc **Zo, nde **No);
#define  Delaunay(Z,N) estiva_Delaunay(&(Z),&(N)) 
extern xyc *estiva_fp2Z(FILE *fp);
#define  fp2Z(fp) estiva_fp2Z(fp) 
extern int main(int argc, char **argv);
#endif
