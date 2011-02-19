#include <stdio.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/solver.h>

typedef long*    I_  ;
typedef long**   I__ ;
typedef double*  D_  ;
typedef double** D__ ;

#define D_  static D_
#define D__ static D__
#define I_  static I_
#define I__ static I__

static int count_NL(void *pA, long N)
{
  long i, j, k, NL ;

  for(NL=0, i=1; i<=N; i++){
    for(k=0, j=1; j<=N; j++) if(mx(pA,i,j) != 0.0) k++;
    if(NL<k) NL = k;
  }
  return NL;
}

static int pcgs();

int estiva_pcgssolver(void* pA, double* x, double* b)
{

  /* (i)  �����η��ȼ���                                                */
  D_     D, B         ;/*    D  �켡������ D(N), B(N)                   */
  D__    A            ;/*    D  �󼡸����� A(N1, 2 * NL)                */
  I__    IA           ;/*    I  �󼡸����� IA(N1, 2 * NL)               */
  D_     R            ;/*    D  �켡������ R(N)                         */
  long   NL, N1, N, ITR, IER                                            ;
  double EPS, S                                                         ;
  D_     X, DD, P, Q  ;/*    D  �켡�������, ���Ǥ� 0��N               */
  I_     M            ;/*    I  �켡������  M(2 * N)  �����            */

  /*--    D����R�ޤǤ�, EPS����S�ޤǤΰ�����, ���֥롼����PCG��Ʊ��   --*/  

  D_     R0, E, H     ;/*    D  �켡������  ���ǿ���N                   */
  D_     W            ;/*    D  �켡������  W(0:N)                      */
  
  long   i, j, k      ;
  
  N   = dim1(b)       ;
  NL  = count_NL(pA,N);
  
  ary1(  D, N   )         ;
  ary2(  A, 2*NL, N+2*NL );
  ary2( IA, 2*NL, N+2*NL );
  ary1(  R, N+1 )         ;
  ary1(  X, N+1 )         ;  ary1( DD, N+1 );  ary1( P, N+1 );  ary1( Q, N+1 ); 
  ary1(  M, 2*N )         ;

  ary1( R0, N   )         ;  ary1(  E, N   );  ary1( H, N   );
  ary1(  W, N+1 )         ;
  

  /* (ii) ��ץ���� �� ���֥롼����                                  */
  /*      ���֥롼�����CALL����Ȥ��ˤ�, �Ĥ����ͤ�Ϳ����.             */
  /*  D   : ����D����1����n���֤˹���A���г����Ǥ�����Ƥ���            */
  /*  N   : ����A�ιԿ�������Ƥ���.                                    */
  /*  N1  : ����A�ιԿ�������Ƥ���. N1��N+2*NL �Ǥʤ��Ȥ����ʤ�.       */
  /*  NL  : ����A�γƹԤˤ������󥼥����ǿ��κ����ͤ�����Ƥ���.        */
  /*  B   : ϢΩ�켡�������α��դ�����Ƥ���.                           */
  /* EPS  : ��«Ƚ���֤�����Ƥ���. �դĤ��� 1.��10^(-7)                */
  /* ITR  : ���ڤ�ޤǤκ��緫�֤����������Ƥ���.                     */

  /*--    D, N, N1, NL, B, EPS, ITR �ˤĤ��Ƥ�, ���֥롼����PCG��Ʊ�� --*/

  /*  A   : ����A�γƹ�1��NL���Ǥ�, ����A�β�������ʬ�γƹԤ��󥼥�     */
  /*        ���Ǥ������. �ޤ�, �ƹԤ�NL+1��2*NL���Ǥ�, �廰����ʬ      */
  /*        �γƹԤ��󥼥����Ǥ������. ������, �ƹԤ��г����Ǥ�����D�� */
  /*        �����.                                                     */
  /* IA   : ����A�����줿���Ǥ����ֹ��, �б�������֤�����Ƥ���.      */
  /*  S   : LUCGSˡ�ΤȤ� 0., MLUCGˡ�ΤȤ���(>0)������Ƥ���.          */  

  for(i=1; i<=N; i++) D[i-1] = mx(pA,i,i); 

  for(i=1; i<=N; i++)for(k=0, j=1; j<i; j++)if(mx(pA,i,j) != 0.0){
    A[k][i-1]  = mx(pA,i,j)  ;
    IA[k][i-1] = j       ;
    k++                  ;
  }

  for(i=1; i<=N; i++)for(k=NL, j=i+1; j<=N; j++)if(mx(pA,i,j) != 0.0){
    A[k][i-1]  = mx(pA,i,j) ;
    IA[k][i-1] = j       ;
    k++                  ;
  }
      

  N1  =  N+2*NL;

  B   = &b[1]  ;
  EPS = 1.0e-7 ;
  ITR = N      ;
  S   = 0.     ;
  
  /* (a) ��󥯤���ˡ                                                   */
  /* CALL PCGS(D,A,IA,N,N1,NL,B,EPS,ITR,S,X,DD,P,Q,R,R0,E,H,W,M,IER)    */

  pcgs(D,A[0],IA[0],&N,&N1,&NL,B,&EPS,&ITR,
			    &S,X,DD,P,Q,R,R0,E,H,W,M,&IER);

  B   = &x[1]  ;
  for(i=0; i<N; i++) B[i] = X[i+1];
  printf("ITR = %ld\n",ITR);
  return IER;
}


#include <math.h>
#define A(i,j)        a[(i)+(j)*dim1]
#define IA(i,j)       ia[(i)+(j)*dim1]
#define forloopL(i,j) for ((j) = 1; (j) <= m[(i)]; (j)++) 
#define forloopU(i,j) for ((j) = *nl + 1; (j) <= *nl + m[*n + (i)]; (j)++) 
#define forall(i)     for ((i) = 1; (i) <= *n; (i)++)
#define mulAij(x)     A(i,j) * (x)[IA(i,j)]

typedef struct {
  long  *ia, *n, *nl, *m, dim1;
  double *a, *d, *dd, *q, *x, *r, *b, *r0, *p, *e, *c1p;
} ILUtype;

typedef long * longp;
typedef double * doublep;

static void ILU(ILUtype ilu)
{
  long   i, j, k, nn, dim1=ilu.dim1;
  double c, ss, sw;
  
  longp ia=ilu.ia, n=ilu.n, nl=ilu.nl, m=ilu.m;
  doublep a=ilu.a, d=ilu.d, dd=ilu.dd, q=ilu.q, x=ilu.x,
    r=ilu.r, b=ilu.b, r0=ilu.r0, p=ilu.p, e=ilu.e;
  
  dd[1] = 1.0 / d[1];

  for (i = 2; i <= *n; i++) {
    ss = d[i];
    for (k = 1; k <= m[i]; k++) {
      nn = IA(i,k);
      for (j = *nl + 1; j <= *nl + m[nn+*n]; j++) 
	if (IA(nn,j) == i) ss -= A(i,k) * A(nn,j) * dd[nn];
    }
    dd[i] = 1.0 / ss;
  }
  
  forall(i) {
    q[i] = d[i] * x[i];
    forloopL(i,j) q[i] += mulAij(x);
    forloopU(i,j) q[i] += mulAij(x);
  }
  
  forall(i) r[i] = b[i] - q[i];
  
  forall(i) {
    forloopL(i,j) r[i] -= mulAij(r);
    r[i] *= dd[i];
  }
  
  for (i = *n; i >= 1; i--) {
    sw = 0.0;
    forloopU(i,j) sw += mulAij(r);
    r[i] -= dd[i] * sw;
  }
  
  c = 0.0;
  forall(i) {
    r0[i] = p[i] = e[i] = r[i];
    c += r[i] * r[i];
  }
  
  *(ilu.c1p) = c;
}


static void presolve(ILUtype ilu, double *q, double *p)
{
  long i, j, dim1=ilu.dim1;
  double sw;

  longp   nl=ilu.nl, n=ilu.n, m=ilu.m, ia=ilu.ia;
  doublep d=ilu.d, dd=ilu.dd, a=ilu.a;

  forall(i) {
    q[i] = d[i] * p[i];
    forloopL(i,j) q[i] += mulAij(p);
    forloopU(i,j) q[i] += mulAij(p);
  }
  
  forall(i) {
    forloopL(i,j) q[i] -= mulAij(q);
    q[i] = dd[i] * q[i];
  }
  
  for (i = *n; i >= 1; i--) {
    sw = 0.0;
    forloopU(i,j) sw += mulAij(q);
    q[i] -= dd[i] * sw;
  }
}


static int pcgs(double *d, double *a, long *ia, 
	long *n, long *n1, long *nl, double *b, double *eps, 
	long *itr, double *s, double *x, double *dd, 
	double *p, double *q, double *r, double *r0, 
	double *e, double *h, double *w, long *m, long *ier)
{
  long dim1, i, j, k;
  double y, c1=0.0, c2, c3, x1, x2, th, res=10000000.0, beta, alpha;
  ILUtype ilu;

  /* Parameter adjustments */
  m--; h--; e--; r0--; b--; d--;
  ia -= (1+*n1);
  a  -= (1+*n1);
  dim1 = *n1;

  /* Function Body */
  *ier = 0;

  if (*n1 < *n || *s < 0.0) { *ier = 2; return 0; }

  th = 1.0;
  if (*s > 0.0 && *s < 1.0) { th = *s; *s = 1.0; }

  for (i = 1; i <= *n << 1; i++) m[i] = 0;

  forall(i) {
    dd[i] = 0.0;
    for (j = 1; j <= *nl; j++)            if (IA(i,j) != 0) m[i]++;
    for (j = *nl + 1; j <= *nl << 1; j++) if (IA(i,j) != 0) m[i + *n]++;
  }


  dd[0] = x[0] = p[0] = q[0] = r[0] = w[0] = 0.0;
  

/*  Incomplete LU Decomposition */

  ilu.ia=ia, ilu.n=n, ilu.nl=nl, ilu.m=m, ilu.dim1=dim1;
  ilu.a=a, ilu.d=d, ilu.dd=dd, ilu.q=q, ilu.x=x,
    ilu.r=r, ilu.b=b, ilu.r0=r0, ilu.p=p, ilu.e=e, ilu.c1p=&c1;

  ILU(ilu);

/*  Iteraton Phase */

    for (k = 1; k <= *itr; k++) {

      presolve(ilu,q,p);

      c2 = 0.0;
      forall(i) c2 += q[i] * r0[i];
      
      if (c2 == 0.0) { *ier = 3; *itr = k; goto finish; }

      alpha = c1 / c2; c3 = x1 = x2 = 0.0;
	
      forall(i) h[i] = e[i] - alpha * q[i];
      forall(i) w[i] = e[i] + h[i];
	
      presolve(ilu,q,w);
      
      forall(i) {
	y       = x[i];
	r[i]   -= alpha * q[i];
	x[i]   += alpha * w[i];
	c3     += r[i] * r0[i];
	x1     += y * y;
	x2     += (x[i]-y) * (x[i]-y);
      }
      
      if (x1 != 0.0) {
	res = sqrt(x2 / x1);
	if (res <= *eps) {
	  *itr = k;
	  *ier = 0;
	  goto finish;
	}
      }

      if (c1 == 0.0) {
	*ier = 4;
	*itr = k;
	goto finish;
      }

      beta = c3 / c1;
      c1 = c3;

      forall(i) {
	e[i] = r[i] + beta * h[i];
	p[i] = e[i] + beta * (h[i] + beta * p[i]);
      }
    }
    *ier = 1;
    
finish:
    *eps = res;

    if (th != 1.0) *s = th;

    return 0;
}
