#include <stdio.h>
#include <string.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/precond.h>
#include <estiva/solver.h>
#include <estiva/op.h>


static MX *A, *AT;
static double *D;

static int matvec(double *alpha, double *x, double *beta, double *y)
{
  static double *t;
  long i, n = A->m; 

  mulmx(t,A,x);
  for(i=0;i<n;i++) t[i] *= *alpha;
  for(i=0;i<n;i++) t[i] += *beta*y[i];
  for(i=0;i<n;i++) y[i] = t[i];
  return 0;
}


static int matvectrans(double *alpha, double *x, double *beta, double *y){
  static double *t;
  long i, n = AT->m;

  mulmx(t,AT,x);
  for(i=0;i<n;i++) t[i] *= *alpha;
  for(i=0;i<n;i++) t[i] += *beta*y[i];
  for(i=0;i<n;i++) y[i] = t[i];
  return 0;
}


int precondjacobi(MX *A, double *x, double *b)
{
  int i,j,J,k,n;
  static double *xold;

  n = dim1(D);

  ary1(xold,n+1);
    
  for (i=0; i<n; i++) xold[i] = 0.0;

  for(k=0;k<3;k++) {
    for(i=0;i<n;i++) {
      x[i]=b[i];
      for(j=0; j< A->n; j++) {
	J = A->IA[i][j];
	if (J != 0) if ( J-1 != i) if(A->A[i][j] !=0.0) {
	  x[i] -= A->A[i][j]*xold[J-1];
	}
      }
      x[i]=x[i]*D[i];
    }
    for(i=0;i<n;i++) xold[i]=x[i];
  }
  return 0;
}


static int psolve(double *x, double *b)
{
  if (!strcmp(getop("-precond"),"none")) { 
    precondnone(A->m,x,b); 
    return 0;
  }
  if (!strcmp(getop("-precond"),"scaling")) { 
    precondscaling(x,D,b); 
    return 0;
  }
  precondjacobi(A,x,b);
  return 0;
}


static int psolvetrans(double *x, double *b)
{  
  if (!strcmp(getop("-precond"),"none")) { 
    precondnone(A->m,x,b); 
    return 0;
  }
  if (!strcmp(getop("-precond"),"scaling")) { 
    precondscaling(x,D,b); 
    return 0;
  }
  precondjacobi(AT,x,b);
  return 0;
}

static int bicg_();

int estiva_bicgsolver(void* pA, double* x, double* b)
{
   int n,ldw,iter,info, i; 
   static double *work, resid;
   printf("phase 1\n");
   A = pA;
   transmx(AT,A);

   printf("phase 2\n");
   n = dim1(b);

   ary1(D,n+1);
   for(i=0;i<n;i++) 
     if (mx(A,i+1,i+1) == 0.0) D[i] = 1.0;
     else                      D[i] = 1.0/mx(A,i+1,i+1);

   ary1(work,n*6+1);
   ldw = n;
   iter = 100 * n;
   resid = 0.0000001;

   for ( i=0; i<n; i++ ) x[i] = b[i];

   printf("phase 3\n");
   bicg_(&n, &b[1], &x[1], &work[1],
	 &ldw, &iter, &resid, matvec, matvectrans, psolve,
	 psolvetrans, &info);
   printf("phase 4\n");
   printf("iter = %d\n",iter);

   return iter;
}

#if 0
int main()
{
  int n;
  static double *b, *x;

  n = 2;
  initmx(A,n+1,n+1); ary1(x,n+1); ary1(b,n+1); 

  mx(A,1,1) = 2.0; mx(A,1,2) = -1.0; b[1] = 10.0;
  mx(A,2,1) = -1.0; mx(A,2,2) = 2.0; b[2] = 4.0;

  estiva_bicgsolver(A,x,b);

  printf("%f\n", x[1] );
  printf("%f\n", x[2] );

  return 0;
}
#endif


/******************* f2c.h ***************************************************/
/* f2c.h  --  Standard Fortran to C header file */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

	- From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

typedef long int integer;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

extern int ungetc();	/* This is not declared in Sun's <stdio.h> */

/* I/O stuff */

#ifdef f2c_i2
/* for -i2 */
typedef short flag;
typedef short ftnlen;
typedef short ftnint;
#else
typedef long flag;
typedef long ftnlen;
typedef long ftnint;
#endif

/*external read, write*/
typedef struct
{	flag cierr;
	ftnint ciunit;
	flag ciend;
	char *cifmt;
	ftnint cirec;
} cilist;

/*internal read, write*/
typedef struct
{	flag icierr;
	char *iciunit;
	flag iciend;
	char *icifmt;
	ftnint icirlen;
	ftnint icirnum;
} icilist;

/*open*/
typedef struct
{	flag oerr;
	ftnint ounit;
	char *ofnm;
	ftnlen ofnmlen;
	char *osta;
	char *oacc;
	char *ofm;
	ftnint orl;
	char *oblnk;
} olist;

/*close*/
typedef struct
{	flag cerr;
	ftnint cunit;
	char *csta;
} cllist;

/*rewind, backspace, endfile*/
typedef struct
{	flag aerr;
	ftnint aunit;
} alist;

/* inquire */
typedef struct
{	flag inerr;
	ftnint inunit;
	char *infile;
	ftnlen infilen;
	ftnint	*inex;	/*parameters in standard's order*/
	ftnint	*inopen;
	ftnint	*innum;
	ftnint	*innamed;
	char	*inname;
	ftnlen	innamlen;
	char	*inacc;
	ftnlen	inacclen;
	char	*inseq;
	ftnlen	inseqlen;
	char 	*indir;
	ftnlen	indirlen;
	char	*infmt;
	ftnlen	infmtlen;
	char	*inform;
	ftnint	informlen;
	char	*inunf;
	ftnlen	inunflen;
	ftnint	*inrecl;
	ftnint	*innrec;
	char	*inblank;
	ftnlen	inblanklen;
} inlist;

#define VOID void

union Multitype {	/* for multiple entry points */
	shortint h;
	integer i;
	real r;
	doublereal d;
	complex c;
	doublecomplex z;
	};

typedef union Multitype Multitype;

typedef long Long;	/* No longer used; formerly in Namelist */

struct Vardesc {	/* for Namelist */
	char *name;
	char *addr;
	ftnlen *dims;
	int  type;
	};
typedef struct Vardesc Vardesc;

struct Namelist {
	char *name;
	Vardesc **vars;
	int nvars;
	};
typedef struct Namelist Namelist;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef shortint (*J_fp)(...);
typedef integer (*I_fp)(...);
typedef real (*R_fp)(...);
typedef doublereal (*D_fp)(...), (*E_fp)(...);
typedef /* Complex */ VOID (*C_fp)(...);
typedef /* Double Complex */ VOID (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
typedef /* Character */ VOID (*H_fp)(...);
typedef /* Subroutine */ int (*S_fp)(...);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef shortint (*J_fp)();
typedef integer (*I_fp)();
typedef real (*R_fp)();
typedef doublereal (*D_fp)(), (*E_fp)();
typedef /* Complex */ VOID (*C_fp)();
typedef /* Double Complex */ VOID (*Z_fp)();
typedef logical (*L_fp)();
typedef shortlogical (*K_fp)();
typedef /* Character */ VOID (*H_fp)();
typedef /* Subroutine */ int (*S_fp)();
#endif
/* E_fp is for real functions when -R is not specified */
typedef VOID C_f;	/* complex function */
typedef VOID H_f;	/* character function */
typedef VOID Z_f;	/* double complex function */
typedef doublereal E_f;	/* real function with -R not specified */

/* undef any lower-case symbols that your C compiler predefines, e.g.: */

#ifndef Skip_f2c_Undefs
#undef cray
#undef gcos
#undef mc68010
#undef mc68020
#undef mips
#undef pdp11
#undef sgi
#undef sparc
#undef sun
#undef sun2
#undef sun3
#undef sun4
#undef u370
#undef u3b
#undef u3b2
#undef u3b5
#undef unix
#undef vax
#endif
#endif


/****************************** static **************************************/


static int daxpy_();
static int dcopy_();
static doublereal ddot_();
static int dlamc1_();
static int dlamc2_();
static doublereal dlamc3_();
static int dlamc4_();
static int dlamc5_();
static doublereal dlamch_();
static doublereal dnrm2_();
static integer do_fio();
static int  e_d();
static integer  e_rsfe();
static integer  e_wsfe();
static int  en_fio();
static int  err__fl();
static int  f__canseek();
static void  f__fatal();
static char*  f__icvt();
static long  f__inode();
static int  f__isdev();
static int  f__nowwriting();
static integer  f_clos();
static void  f_exit();
static void  f_init();
static char*  f_list();
static integer  f_open();
static char*  f_s();
static int  fk_open();
static void  fmt_bg();
static doublereal  getbreak_();
static char*  gt_num();
static char*  i_tem();
static logical  lsame_();
static int  mv_cur();
static int  ne_d();
static int  op_gen();
static int  pars_f();
static doublereal  pow_di();
static integer  s_wsfe();
static void  sig_die();
static int  t_runc();
static int  type_f();
static int  w_ed();
static int  w_ned();
static int  wrt_E();
static int  wrt_F();
static int  wrt_L();
static int  x_putc();
static int  x_wSL();
static int  xw_end();
static int  xw_rev();


#ifdef NON_ANSI_RW_MODES
char *f__r_mode[2] = {"r", "r"};
static char *f__w_mode[4] = {"w", "w", "r+w", "r+w"};
#else
static char *f__r_mode[2] = {"rb", "r"};
static char *f__w_mode[4] = {"wb", "w", "r+b", "r+"};
#endif


/*error messages*/
static char *F_err[] =
{
	"error in format",				/* 100 */
	"illegal unit number",				/* 101 */
	"formatted io not allowed",			/* 102 */
	"unformatted io not allowed",			/* 103 */
	"direct io not allowed",			/* 104 */
	"sequential io not allowed",			/* 105 */
	"can't backspace file",				/* 106 */
	"null file name",				/* 107 */
	"can't stat file",				/* 108 */
	"unit not connected",				/* 109 */
	"off end of record",				/* 110 */
	"truncation failed in endfile",			/* 111 */
	"incomprehensible list input",			/* 112 */
	"out of free space",				/* 113 */
	"unit not connected",				/* 114 */
	"read unexpected character",			/* 115 */
	"bad logical input field",			/* 116 */
	"bad variable type",				/* 117 */
	"bad namelist name",				/* 118 */
	"variable not in namelist",			/* 119 */
	"no end record",				/* 120 */
	"variable count incorrect",			/* 121 */
	"subscript for scalar variable",		/* 122 */
	"invalid array section",			/* 123 */
	"substring out of bounds",			/* 124 */
	"subscript out of bounds",			/* 125 */
	"can't read file",				/* 126 */
	"can't write file",				/* 127 */
	"'new' file exists",				/* 128 */
	"can't append to file"				/* 129 */
};


/********************************** fio.h ***********************************/
#include "stdio.h"
#include "errno.h"
#ifndef NULL
/* ANSI C */
#include "stddef.h"
#endif

#ifndef SEEK_SET
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

#ifdef MSDOS
#ifndef NON_UNIX_STDIO
#define NON_UNIX_STDIO
#endif
#endif

#ifdef UIOLEN_int
typedef int uiolen;
#else
typedef long uiolen;
#endif

/*units*/
typedef struct
{	FILE *ufd;	/*0=unconnected*/
	char *ufnm;
#ifndef MSDOS
	long uinode;
	int udev;
#endif
	int url;	/*0=sequential*/
	flag useek;	/*true=can backspace, use dir, ...*/
	flag ufmt;
	flag uprnt;
	flag ublnk;
	flag uend;
	flag uwrt;	/*last io was write*/
	flag uscrtch;
} unit;

static flag f__init;
static cilist *f__elist;	/*active external io list*/
static flag f__reading,f__external,f__sequential,f__formatted;
#undef Void
#ifdef KR_headers
#define Void /*void*/
static int (*f__getn)(),(*f__putn)();	/*for formatted io*/
extern long f__inode();
extern VOID sig_die();
extern int (*f__donewrec)(), t_putc(), x_wSL();
static int c_sfe(), err_fl(), xrd_SL();
#else
#define Void void
#ifdef __cplusplus
extern "C" {
#endif
static int (*f__putn)(int);	/*for formatted io*/
extern long f__inode(char*,int*);
extern void sig_die(char*,int);
extern void f__fatal(int,char*);
extern int t_runc(alist*);
extern int  f__nowwriting(unit*);
extern int fk_open(int,int,ftnint);
extern int en_fio(void);
extern void f_init(void);
static int (*f__donewrec)(void), x_wSL(void);
static void g_char(char*,ftnlen,char*);
static int c_sfe(cilist*);
extern int isatty(int);
extern int err__fl(int,int,char*);
extern int xrd_SL(void);
#ifdef __cplusplus
	}
#endif
#endif
static int (*f__doend)(Void);
static FILE *f__cf;	/*current file*/
static unit *f__curunit;	/*current unit*/
static unit f__units[];
#define err(f,m,s) {if(f) errno= m; else f__fatal(m,s); return(m);}
#define errfl(f,m,s) return err__fl(f,m,s)

/*Table sizes*/
#define MXUNIT 100

static int f__recpos;	/*position in current record*/
static int f__cursor;	/* offset to move to */
static int f__hiwater;	/* so TL doesn't confuse us */

#define WRITE	1
#define READ	2
#define SEQ	3
#define DIR	4
#define FMT	5
#define UNF	6
#define EXT	7
#define INT	8

#define buf_end(x) (x->_flag & _IONBF ? x->_ptr : x->_base + BUFSIZ)


/******************************** fmt.h **************************************/
struct f__syl
{	int op,p1,p2,p3;
};
#define RET1 1
#define REVERT 2
#define GOTO 3
#define X 4
#define SLASH 5
#define STACK 6
#define I 7
#define ED 8
#define NED 9
#define IM 10
#define APOS 11
#define H 12
#define TL 13
#define TR 14
#define T 15
#define COLON 16
#define S 17
#define SP 18
#define SS 19
#define P 20
#define BN 21
#define BZ 22
#define F 23
#define E 24
#define EE 25
#define D 26
#define G 27
#define GE 28
#define L 29
#define A 30
#define AW 31
#define O 32
#define NONL 33
#define OM 34
#define Z 35
#define ZM 36
static struct f__syl f__syl[];
static int f__pc,f__parenlvl,f__revloc;
typedef union
{	real pf;
	doublereal pd;
} ufloat;
typedef union
{	short is;
	char ic;
	long il;
#ifdef Allow_TYQUAD
	longint ili;
#endif
} Uint;
#ifdef KR_headers
static int (*f__doed)(),(*f__doned)();
static int (*f__dorevert)();
extern int rd_ed(),rd_ned();
extern int w_ed(),w_ned();
#else
#ifdef __cplusplus
extern "C" {
#endif
static int (*f__doed)(struct f__syl*, char*, ftnlen),(*f__doned)(struct f__syl*);
static int (*f__dorevert)(void);
extern void fmt_bg(void);
extern int pars_f(char*);
extern int rd_ed(struct f__syl*, char*, ftnlen),rd_ned(struct f__syl*);
extern int w_ed(struct f__syl*, char*, ftnlen),w_ned(struct f__syl*);
extern int wrt_E(ufloat*, int, int, int, ftnlen);
extern int wrt_F(ufloat*, int, int, ftnlen);
extern int wrt_L(Uint*, int, ftnlen);
#ifdef __cplusplus
	}
#endif
#endif
static flag f__cblank,f__cplus,f__workdone, f__nonl;
static char *f__fmtbuf;
static int f__scale;
#define GET(x) if((x=(*f__getn)())<0) return(x)
#define VAL(x) (x!='\n'?x:' ')
#define PUT(x) (*f__putn)(x)
extern int f__cursor;

/******************************** rawid.h ************************************/
#ifdef KR_headers
extern FILE *fdopen();
#else
#ifdef MSDOS
#include "io.h"
#define close _close
#define creat _creat
#define open _open
#define read _read
#define write _write
#endif
#ifdef __cplusplus
extern "C" {
#endif
#ifndef MSDOS
#ifdef OPEN_DECL
extern int creat(const char*,int), open(const char*,int);
#endif
extern int close(int);
extern int read(int,void*,size_t), write(int,void*,size_t);
extern int unlink(const char*);
#ifndef _POSIX_SOURCE
#ifndef NON_UNIX_STDIO
/* 
 * This type declaration may not be consistent with the one in <stdio.h>.
 * extern FILE *fdopen(int, const char*); 
 */
//extern FILE *fdopen(int, char*);
#endif
#endif
#endif

extern char *mktemp(char*);

#ifdef __cplusplus
	}
#endif
#endif

#include "fcntl.h"

#ifndef O_WRONLY
#define O_RDONLY 0
#define O_WRONLY 1
#endif

/********************************* fp.h **************************************/
#define FMAX 40
#define EXPMAXDIGS 8
#define EXPMAX 99999999
/* FMAX = max number of nonzero digits passed to atof() */
/* EXPMAX = 10^EXPMAXDIGS - 1 = largest allowed exponent absolute value */



/* MAXFRACDIGS and MAXINTDIGS are for wrt_F -- bounds (not necessarily
   tight) on the maximum number of digits to the right and left of
 * the decimal point.
 */

#ifdef VAX
#define MAXFRACDIGS 56
#define MAXINTDIGS 38
#else
#ifdef CRAY
#define MAXFRACDIGS 9880
#define MAXINTDIGS 9864
#else
/* values that suffice for IEEE double */
#define MAXFRACDIGS 344
#define MAXINTDIGS 308
#endif
#endif

/*****************************************************************************/

doublereal getbreak_()
{
  /* System generated locals */
  doublereal ret_val, d__1;

  /* Local variables */
  extern doublereal dlamch_();
  static doublereal eps;


  /*     Get breakdown parameter tolerance; for the test routine, */
  /*     set to machine precision. */


  eps = dlamch_("EPS", 3L);
  /* Computing 2nd power */
  d__1 = eps;
  ret_val = d__1 * d__1;

  return ret_val;

} /* getbreak_ */


/********************************* BiCG.c ************************************/
/* BiCG.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Table of constant values */

static integer c__BiCG1 = 1;
static doublereal c_b5 = -1.;
static doublereal c_b6 = 0.;
static doublereal c_b28 = 1.;

/*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the
*     Solution of Linear Systems: Building Blocks for Iterative
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*  Purpose
*  =======
*
*  BiCG solves the linear system Ax = b using the
*  BiConjugate Gradient iterative method with preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER.
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,6).
*          Workspace for residual, direction vector, etc.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  MATVECTRANS  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A'*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A' is the tranpose of a matrix A. Vector x must remain
*          unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVECTRANS( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must
*          remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*         The preconditioner is passed into the routine in a common block
*
*  PSOLVETRANS  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M'*x = b,
*
*          where x and y are vectors, and M' is the tranpose of a
*          matrix M. Vector b must remain unchanged. The solution 
*          is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVETRANS( X, B )
*
*         The preconditioner is passed into the routine in a common block
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter.
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO < BREAKTOL: RHO and RTLD have become
*                                       orthogonal.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  BLAS CALLS:   DAXPY, DCOPY, DDOT, DNRM2
*  ============================================================ */

static int bicg_(n, b, x, work, ldw, iter, resid, matvec, matvectrans, 
          psolve, psolvetrans, info)
   integer *n, *ldw, *iter, *info;
   doublereal *b, *x, *work, *resid;
   int (*matvec) (), (*matvectrans) (), (*psolve) (), (*psolvetrans) ();
{
    /* System generated locals */
    integer work_dim1, work_offset;
    doublereal d__1;

    /* Local variables */
    static doublereal beta;
    extern doublereal ddot_();
    static integer ptld, rtld, qtld;
    extern doublereal getbreak_();
    static integer ztld;
    static doublereal bnrm2;
    extern doublereal dnrm2_();
    static integer p, q, r;
    static doublereal alpha;
    static integer z;
    extern /* Subroutine */ int dcopy_();
    static integer maxit;
    extern /* Subroutine */ int daxpy_();
    static doublereal rhotol, rho, tol, rho1;

    /* Parameter adjustments */
    work_dim1 = *ldw;
    work_offset = work_dim1 + 1;
    work -= work_offset;
    --x;
    --b;

/*  Executable Statements */

    *info = 0;

/*     Test the input parameters. */

    if (*n < 0) {
	*info = -1;
    } else if (*ldw < max(1,*n)) {
	*info = -2;
    } else if (*iter <= 0) {
	*info = -3;
    }
    if (*info != 0) {
	return 0;
    }

    maxit = *iter;
    tol = *resid;

/*     Alias workspace columns. */

    r = 1;
    rtld = 2;
    z = 3;
    ztld = 4;
    p = 5;
    ptld = 6;
    q = 3;
    qtld = 4;

/*     Set breakdown parameters. */

    rhotol = getbreak_();

/*     Set initial residual. */

    dcopy_(n, &b[1], &c__BiCG1, &work[r * work_dim1 + 1], &c__BiCG1);
    if (dnrm2_(n, &x[1], &c__BiCG1) != 0.) {
	(*matvec)(&c_b5, &x[1], &c_b28, &work[r * work_dim1 + 1]);
	if (dnrm2_(n, &work[r * work_dim1 + 1], &c__BiCG1) <= tol) {
	    goto L30;
	}
    }
    dcopy_(n, &work[r * work_dim1 + 1], &c__BiCG1, &work[rtld * work_dim1 + 1], &
	    c__BiCG1);
    bnrm2 = dnrm2_(n, &b[1], &c__BiCG1);
    if (bnrm2 == 0.) {
	bnrm2 = 1.;
    }

    *iter = 0;

L10:

/*     Perform BiConjugate Gradient iteration. */

    ++(*iter);

/*        Compute direction vectors PK and PTLD. */

    (*psolve)(&work[z * work_dim1 + 1], &work[r * work_dim1 + 1]);
    (*psolvetrans)(&work[ztld * work_dim1 + 1], &work[rtld * work_dim1 + 1]);

    rho = ddot_(n, &work[z * work_dim1 + 1], &c__BiCG1, &work[rtld * work_dim1 + 
	    1], &c__BiCG1);
    if (abs(rho) < rhotol) {
	goto L25;
    }

    if (*iter > 1) {
	beta = rho / rho1;
	daxpy_(n, &beta, &work[p * work_dim1 + 1], &c__BiCG1, &work[z * work_dim1 
		+ 1], &c__BiCG1);
	daxpy_(n, &beta, &work[ptld * work_dim1 + 1], &c__BiCG1, &work[ztld * 
		work_dim1 + 1], &c__BiCG1);
	dcopy_(n, &work[z * work_dim1 + 1], &c__BiCG1, &work[p * work_dim1 + 1], &
		c__BiCG1);
	dcopy_(n, &work[ztld * work_dim1 + 1], &c__BiCG1, &work[ptld * work_dim1 
		+ 1], &c__BiCG1);
    } else {
	dcopy_(n, &work[z * work_dim1 + 1], &c__BiCG1, &work[p * work_dim1 + 1], &
		c__BiCG1);
	dcopy_(n, &work[ztld * work_dim1 + 1], &c__BiCG1, &work[ptld * work_dim1 
		+ 1], &c__BiCG1);
    }

    (*matvec)(&c_b28, &work[p * work_dim1 + 1], &c_b6, &work[q * work_dim1 + 
	    1]);
    (*matvectrans)(&c_b28, &work[ptld * work_dim1 + 1], &c_b6, &work[qtld * 
	    work_dim1 + 1]);
    alpha = rho / ddot_(n, &work[ptld * work_dim1 + 1], &c__BiCG1, &work[q * 
	    work_dim1 + 1], &c__BiCG1);

/*        Compute current solution vector x. */

    daxpy_(n, &alpha, &work[p * work_dim1 + 1], &c__BiCG1, &x[1], &c__BiCG1);

/*        Compute residual vector rk, find norm, */
/*        then check for tolerance. */

    d__1 = -alpha;
    daxpy_(n, &d__1, &work[q * work_dim1 + 1], &c__BiCG1, &work[r * work_dim1 + 1]
	    , &c__BiCG1);

    *resid = dnrm2_(n, &work[r * work_dim1 + 1], &c__BiCG1) / bnrm2;
    if (*resid <= tol) {
	goto L30;
    }
    if (*iter == maxit) {
	goto L20;
    }

    d__1 = -alpha;
    daxpy_(n, &d__1, &work[qtld * work_dim1 + 1], &c__BiCG1, &work[rtld * 
	    work_dim1 + 1], &c__BiCG1);
    rho1 = rho;

    goto L10;

L20:

/*     Iteration fails. */

    *info = 1;
    return 0;

L25:

/*     Set breakdown flag. */

    *info = -10;
    return 0;

L30:

/*     Iteration successful; return. */

    return 0;

/*     End of BICG */

}



/********************************* dlapack.c *********************************/

/* dlapack.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* Table of constant values */

static integer c__1 = 1; /****************************************************/
static doublereal c_b21 = 0.;









doublereal dlamch_(cmach, cmach_len)
char *cmach;
ftnlen cmach_len;
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double pow_di();

    /* Local variables */
    static doublereal base;
    static integer beta;
    static doublereal emin, prec, emax;
    static integer imin, imax;
    static logical lrnd;
    static doublereal rmin, rmax, t, rmach;
    extern logical lsame_();
    static doublereal small, sfmin;
    extern /* Subroutine */ int dlamc2_();
    static integer it;
    static doublereal rnd, eps;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMCH determines double precision machine parameters. */

/*  Arguments */
/*  ========= */

/*  CMACH   (input) CHARACTER*1 */
/*          Specifies the value to be returned by DLAMCH: */
/*          = 'E' or 'e',   DLAMCH := eps */
/*          = 'S' or 's ,   DLAMCH := sfmin */
/*          = 'B' or 'b',   DLAMCH := base */
/*          = 'P' or 'p',   DLAMCH := eps*base */
/*          = 'N' or 'n',   DLAMCH := t */
/*          = 'R' or 'r',   DLAMCH := rnd */
/*          = 'M' or 'm',   DLAMCH := emin */
/*          = 'U' or 'u',   DLAMCH := rmin */
/*          = 'L' or 'l',   DLAMCH := emax */
/*          = 'O' or 'o',   DLAMCH := rmax */

/*          where */

/*          eps   = relative machine precision */
/*          sfmin = safe minimum, such that 1/sfmin does not overflow */
/*          base  = base of the machine */
/*          prec  = eps*base */
/*          t     = number of (base) digits in the mantissa */
/*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise */
/*          emin  = minimum exponent before (gradual) underflow */
/*          rmin  = underflow threshold - base**(emin-1) */
/*          emax  = largest exponent before overflow */
/*          rmax  = overflow threshold  - (base**emax)*(1-eps) */

/* ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */

    if (first) {
	first = FALSE_;
	dlamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
	base = (doublereal) beta;
	t = (doublereal) it;
	if (lrnd) {
	    rnd = 1.;
	    i__1 = 1 - it;
	    eps = pow_di(&base, &i__1) / 2;
	} else {
	    rnd = 0.;
	    i__1 = 1 - it;
	    eps = pow_di(&base, &i__1);
	}
	prec = eps * base;
	emin = (doublereal) imin;
	emax = (doublereal) imax;
	sfmin = rmin;
	small = 1. / rmax;
	if (small >= sfmin) {

/*           Use SMALL plus a bit, to avoid the possibility of rou
nding */
/*           causing overflow when computing  1/sfmin. */

	    sfmin = small * (eps + 1.);
	}
    }

    if (lsame_(cmach, "E", 1L, 1L)) {
	rmach = eps;
    } else if (lsame_(cmach, "S", 1L, 1L)) {
	rmach = sfmin;
    } else if (lsame_(cmach, "B", 1L, 1L)) {
	rmach = base;
    } else if (lsame_(cmach, "P", 1L, 1L)) {
	rmach = prec;
    } else if (lsame_(cmach, "N", 1L, 1L)) {
	rmach = t;
    } else if (lsame_(cmach, "R", 1L, 1L)) {
	rmach = rnd;
    } else if (lsame_(cmach, "M", 1L, 1L)) {
	rmach = emin;
    } else if (lsame_(cmach, "U", 1L, 1L)) {
	rmach = rmin;
    } else if (lsame_(cmach, "L", 1L, 1L)) {
	rmach = emax;
    } else if (lsame_(cmach, "O", 1L, 1L)) {
	rmach = rmax;
    }

    ret_val = rmach;
    return ret_val;

/*     End of DLAMCH */

} /* dlamch_ */


/* *********************************************************************** */

/* Subroutine */ int dlamc1_(beta, t, rnd, ieee1)
integer *beta, *t;
logical *rnd, *ieee1;
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static logical lrnd;
    static doublereal a, b, c, f;
    static integer lbeta;
    static doublereal savec;
    extern doublereal dlamc3_();
    static logical lieee1;
    static doublereal t1, t2;
    static integer lt;
    static doublereal one, qtr;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC1 determines the machine parameters given by BETA, T, RND, and */
/*  IEEE1. */

/*  Arguments */
/*  ========= */

/*  BETA    (output) INTEGER */
/*          The base of the machine. */

/*  T       (output) INTEGER */
/*          The number of ( BETA ) digits in the mantissa. */

/*  RND     (output) LOGICAL */
/*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
/*          chopping  ( RND = .FALSE. )  occurs in addition. This may not 
*/
/*          be a reliable guide to the way in which the machine performs 
*/
/*          its arithmetic. */

/*  IEEE1   (output) LOGICAL */
/*          Specifies whether rounding appears to be done in the IEEE */
/*          'round to nearest' style. */

/*  Further Details */
/*  =============== */

/*  The routine is based on the routine  ENVRON  by Malcolm and */
/*  incorporates suggestions by Gentleman and Marovich. See */

/*     Malcolm M. A. (1972) Algorithms to reveal properties of */
/*        floating-point arithmetic. Comms. of the ACM, 15, 949-951. */

/*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms */
/*        that reveal properties of floating point arithmetic units. */
/*        Comms. of the ACM, 17, 276-277. */

/* ===================================================================== 
*/

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */

    if (first) {
	first = FALSE_;
	one = 1.;

/*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BE
TA, */
/*        IEEE1, T and RND. */

/*        Throughout this routine  we use the function  DLAMC3  to ens
ure */
/*        that relevant values are  stored and not held in registers, 
 or */
/*        are not affected by optimizers. */

/*        Compute  a = 2.0**m  with the  smallest positive integer m s
uch */
/*        that */

/*           fl( a + 1.0 ) = a. */

	a = 1.;
	c = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
L10:
	if (c == one) {
	    a *= 2;
	    c = dlamc3_(&a, &one);
	    d__1 = -a;
	    c = dlamc3_(&c, &d__1);
	    goto L10;
	}
/* +       END WHILE */

/*        Now compute  b = 2.0**m  with the smallest positive integer 
m */
/*        such that */

/*           fl( a + b ) .gt. a. */

	b = 1.;
	c = dlamc3_(&a, &b);

/* +       WHILE( C.EQ.A )LOOP */
L20:
	if (c == a) {
	    b *= 2;
	    c = dlamc3_(&a, &b);
	    goto L20;
	}
/* +       END WHILE */

/*        Now compute the base.  a and c  are neighbouring floating po
int */
/*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and
 so */
/*        their difference is beta. Adding 0.25 to c is to ensure that
 it */
/*        is truncated to beta and not ( beta - 1 ). */

	qtr = one / 4;
	savec = c;
	d__1 = -a;
	c = dlamc3_(&c, &d__1);
	lbeta = (integer) (c + qtr);

/*        Now determine whether rounding or chopping occurs,  by addin
g a */
/*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to 
 a. */

	b = (doublereal) lbeta;
	d__1 = b / 2;
	d__2 = -b / 100;
	f = dlamc3_(&d__1, &d__2);
	c = dlamc3_(&f, &a);
	if (c == a) {
	    lrnd = TRUE_;
	} else {
	    lrnd = FALSE_;
	}
	d__1 = b / 2;
	d__2 = b / 100;
	f = dlamc3_(&d__1, &d__2);
	c = dlamc3_(&f, &a);
	if (lrnd && c == a) {
	    lrnd = FALSE_;
	}

/*        Try and decide whether rounding is done in the  IEEE  'round
 to */
/*        nearest' style. B/2 is half a unit in the last place of the 
two */
/*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  
bit */
/*        zero, and SAVEC is odd. Thus adding B/2 to A should not  cha
nge */
/*        A, but adding B/2 to SAVEC should change SAVEC. */

	d__1 = b / 2;
	t1 = dlamc3_(&d__1, &a);
	d__1 = b / 2;
	t2 = dlamc3_(&d__1, &savec);
	lieee1 = t1 == a && t2 > savec && lrnd;

/*        Now find  the  mantissa, t.  It should  be the  integer part
 of */
/*        log to the base beta of a,  however it is safer to determine
  t */
/*        by powering.  So we find t as the smallest positive integer 
for */
/*        which */

/*           fl( beta**t + 1.0 ) = 1.0. */

	lt = 0;
	a = 1.;
	c = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
L30:
	if (c == one) {
	    ++lt;
	    a *= lbeta;
	    c = dlamc3_(&a, &one);
	    d__1 = -a;
	    c = dlamc3_(&c, &d__1);
	    goto L30;
	}
/* +       END WHILE */

    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *ieee1 = lieee1;
    return 0;

/*     End of DLAMC1 */

} /* dlamc1_ */


/* *********************************************************************** */

/* Subroutine */ int dlamc2_(beta, t, rnd, eps, emin, rmin, emax, rmax)
integer *beta, *t;
logical *rnd;
doublereal *eps;
integer *emin;
doublereal *rmin;
integer *emax;
doublereal *rmax;
{
    /* Initialized data */

    static logical first = TRUE_;
    static logical iwarn = FALSE_;

    /* Format strings */
    static char fmt_9999[] = "(//\002 WARNING. The value EMIN may be incorre\
ct:-\002,\002  EMIN = \002,i8,/\002 If, after inspection, the value EMIN loo\
ks\002,\002 acceptable please comment out \002,/\002 the IF block as marked \
within the code of routine\002,\002 DLAMC2,\002,/\002 otherwise supply EMIN \
explicitly.\002,/)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    double pow_di();
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static logical ieee;
    static doublereal half;
    static logical lrnd;
    static doublereal leps, zero, a, b, c;
    static integer i, lbeta;
    static doublereal rbase;
    static integer lemin, lemax, gnmin;
    static doublereal small;
    static integer gpmin;
    static doublereal third, lrmin, lrmax, sixth;
    extern /* Subroutine */ int dlamc1_();
    extern doublereal dlamc3_();
    static logical lieee1;
    extern /* Subroutine */ int dlamc4_(), dlamc5_();
    static integer lt, ngnmin, ngpmin;
    static doublereal one, two;

    /* Fortran I/O blocks */
    static cilist io___231 = { 0, 6, 0, fmt_9999, 0 };



/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC2 determines the machine parameters specified in its argument */
/*  list. */

/*  Arguments */
/*  ========= */

/*  BETA    (output) INTEGER */
/*          The base of the machine. */

/*  T       (output) INTEGER */
/*          The number of ( BETA ) digits in the mantissa. */

/*  RND     (output) LOGICAL */
/*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
/*          chopping  ( RND = .FALSE. )  occurs in addition. This may not 
*/
/*          be a reliable guide to the way in which the machine performs 
*/
/*          its arithmetic. */

/*  EPS     (output) DOUBLE PRECISION */
/*          The smallest positive number such that */

/*             fl( 1.0 - EPS ) .LT. 1.0, */

/*          where fl denotes the computed value. */

/*  EMIN    (output) INTEGER */
/*          The minimum exponent before (gradual) underflow occurs. */

/*  RMIN    (output) DOUBLE PRECISION */
/*          The smallest normalized number for the machine, given by */
/*          BASE**( EMIN - 1 ), where  BASE  is the floating point value 
*/
/*          of BETA. */

/*  EMAX    (output) INTEGER */
/*          The maximum exponent before overflow occurs. */

/*  RMAX    (output) DOUBLE PRECISION */
/*          The largest positive number for the machine, given by */
/*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point 
*/
/*          value of BETA. */

/*  Further Details */
/*  =============== */

/*  The computation of  EPS  is based on a routine PARANOIA by */
/*  W. Kahan of the University of California at Berkeley. */

/* ===================================================================== 
*/

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */

    if (first) {
	first = FALSE_;
	zero = 0.;
	one = 1.;
	two = 2.;

/*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values
 of */
/*        BETA, T, RND, EPS, EMIN and RMIN. */

/*        Throughout this routine  we use the function  DLAMC3  to ens
ure */
/*        that relevant values are stored  and not held in registers, 
 or */
/*        are not affected by optimizers. */

/*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. 
*/

	dlamc1_(&lbeta, &lt, &lrnd, &lieee1);

/*        Start to find EPS. */

	b = (doublereal) lbeta;
	i__1 = -lt;
	a = pow_di(&b, &i__1);
	leps = a;

/*        Try some tricks to see whether or not this is the correct  E
PS. */

	b = two / 3;
	half = one / 2;
	d__1 = -half;
	sixth = dlamc3_(&b, &d__1);
	third = dlamc3_(&sixth, &sixth);
	d__1 = -half;
	b = dlamc3_(&third, &d__1);
	b = dlamc3_(&b, &sixth);
	b = abs(b);
	if (b < leps) {
	    b = leps;
	}

	leps = 1.;

/* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
L10:
	if (leps > b && b > zero) {
	    leps = b;
	    d__1 = half * leps;
/* Computing 5th power */
	    d__3 = two, d__4 = d__3, d__3 *= d__3;
/* Computing 2nd power */
	    d__5 = leps;
	    d__2 = d__4 * (d__3 * d__3) * (d__5 * d__5);
	    c = dlamc3_(&d__1, &d__2);
	    d__1 = -c;
	    c = dlamc3_(&half, &d__1);
	    b = dlamc3_(&half, &c);
	    d__1 = -b;
	    c = dlamc3_(&half, &d__1);
	    b = dlamc3_(&half, &c);
	    goto L10;
	}
/* +       END WHILE */

	if (a < leps) {
	    leps = a;
	}

/*        Computation of EPS complete. */

/*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3
)). */
/*        Keep dividing  A by BETA until (gradual) underflow occurs. T
his */
/*        is detected when we cannot recover the previous A. */

	rbase = one / lbeta;
	small = one;
	for (i = 1; i <= 3; ++i) {
	    d__1 = small * rbase;
	    small = dlamc3_(&d__1, &zero);
/* L20: */
	}
	a = dlamc3_(&one, &small);
	dlamc4_(&ngpmin, &one, &lbeta);
	d__1 = -one;
	dlamc4_(&ngnmin, &d__1, &lbeta);
	dlamc4_(&gpmin, &a, &lbeta);
	d__1 = -a;
	dlamc4_(&gnmin, &d__1, &lbeta);
	ieee = FALSE_;

	if (ngpmin == ngnmin && gpmin == gnmin) {
	    if (ngpmin == gpmin) {
		lemin = ngpmin;
/*            ( Non twos-complement machines, no gradual under
flow; */
/*              e.g.,  VAX ) */
	    } else if (gpmin - ngpmin == 3) {
		lemin = ngpmin - 1 + lt;
		ieee = TRUE_;
/*            ( Non twos-complement machines, with gradual und
erflow; */
/*              e.g., IEEE standard followers ) */
	    } else {
		lemin = min(ngpmin,gpmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else if (ngpmin == gpmin && ngnmin == gnmin) {
	    if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1) {
		lemin = max(ngpmin,ngnmin);
/*            ( Twos-complement machines, no gradual underflow
; */
/*              e.g., CYBER 205 ) */
	    } else {
		lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1 && gpmin == gnmin)
		 {
	    if (gpmin - min(ngpmin,ngnmin) == 3) {
		lemin = max(ngpmin,ngnmin) - 1 + lt;
/*            ( Twos-complement machines with gradual underflo
w; */
/*              no known machine ) */
	    } else {
		lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else {
/* Computing MIN */
	    i__1 = min(ngpmin,ngnmin), i__1 = min(i__1,gpmin);
	    lemin = min(i__1,gnmin);
/*         ( A guess; no known machine ) */
	    iwarn = TRUE_;
	}
/* ** */
/* Comment out this if block if EMIN is ok */
	if (iwarn) {
	    first = TRUE_;
	    s_wsfe(&io___231);
	    do_fio(&c__1, (char *)&lemin, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
/* ** */

/*        Assume IEEE arithmetic if we found denormalised  numbers abo
ve, */
/*        or if arithmetic seems to round in the  IEEE style,  determi
ned */
/*        in routine DLAMC1. A true IEEE machine should have both  thi
ngs */
/*        true; however, faulty machines may have one or the other. */

	ieee = ieee || lieee1;

/*        Compute  RMIN by successive division by  BETA. We could comp
ute */
/*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow dur
ing */
/*        this computation. */

	lrmin = 1.;
	i__1 = 1 - lemin;
	for (i = 1; i <= i__1; ++i) {
	    d__1 = lrmin * rbase;
	    lrmin = dlamc3_(&d__1, &zero);
/* L30: */
	}

/*        Finally, call DLAMC5 to compute EMAX and RMAX. */

	dlamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *eps = leps;
    *emin = lemin;
    *rmin = lrmin;
    *emax = lemax;
    *rmax = lrmax;

    return 0;


/*     End of DLAMC2 */

} /* dlamc2_ */


/* *********************************************************************** */

doublereal dlamc3_(a, b)
doublereal *a, *b;
{
    /* System generated locals */
    doublereal ret_val;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing 
*/
/*  the addition of  A  and  B ,  for use in situations where optimizers 
*/
/*  might hold one of these in a register. */

/*  Arguments */
/*  ========= */

/*  A, B    (input) DOUBLE PRECISION */
/*          The values A and B. */

/* ===================================================================== 
*/

/*     .. Executable Statements .. */

    ret_val = *a + *b;

    return ret_val;

/*     End of DLAMC3 */

} /* dlamc3_ */


/* *********************************************************************** */

/* Subroutine */ int dlamc4_(emin, start, base)
integer *emin;
doublereal *start;
integer *base;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal zero, a;
    static integer i;
    static doublereal rbase, b1, b2, c1, c2, d1, d2;
    extern doublereal dlamc3_();
    static doublereal one;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC4 is a service routine for DLAMC2. */

/*  Arguments */
/*  ========= */

/*  EMIN    (output) EMIN */
/*          The minimum exponent before (gradual) underflow, computed by 
*/
/*          setting A = START and dividing by BASE until the previous A */
/*          can not be recovered. */

/*  START   (input) DOUBLE PRECISION */
/*          The starting point for determining EMIN. */

/*  BASE    (input) INTEGER */
/*          The base of the machine. */

/* ===================================================================== 
*/

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    a = *start;
    one = 1.;
    rbase = one / *base;
    zero = 0.;
    *emin = 1;
    d__1 = a * rbase;
    b1 = dlamc3_(&d__1, &zero);
    c1 = a;
    c2 = a;
    d1 = a;
    d2 = a;
/* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND. */
/*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
L10:
    if (c1 == a && c2 == a && d1 == a && d2 == a) {
	--(*emin);
	a = b1;
	d__1 = a / *base;
	b1 = dlamc3_(&d__1, &zero);
	d__1 = b1 * *base;
	c1 = dlamc3_(&d__1, &zero);
	d1 = zero;
	i__1 = *base;
	for (i = 1; i <= i__1; ++i) {
	    d1 += b1;
/* L20: */
	}
	d__1 = a * rbase;
	b2 = dlamc3_(&d__1, &zero);
	d__1 = b2 / rbase;
	c2 = dlamc3_(&d__1, &zero);
	d2 = zero;
	i__1 = *base;
	for (i = 1; i <= i__1; ++i) {
	    d2 += b2;
/* L30: */
	}
	goto L10;
    }
/* +    END WHILE */

    return 0;

/*     End of DLAMC4 */

} /* dlamc4_ */


/* *********************************************************************** */

/* Subroutine */ int dlamc5_(beta, p, emin, ieee, emax, rmax)
integer *beta, *p, *emin;
logical *ieee;
integer *emax;
doublereal *rmax;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer lexp;
    static doublereal oldy;
    static integer uexp, i;
    static doublereal y, z;
    static integer nbits;
    extern doublereal dlamc3_();
    static doublereal recbas;
    static integer exbits, expsum, try_;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC5 attempts to compute RMAX, the largest machine floating-point */
/*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum */
/*  approximately to a power of 2.  It will fail on machines where this */
/*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625, 
*/
/*  EMAX = 28718).  It will also fail if the value supplied for EMIN is */
/*  too large (i.e. too close to zero), probably with overflow. */

/*  Arguments */
/*  ========= */

/*  BETA    (input) INTEGER */
/*          The base of floating-point arithmetic. */

/*  P       (input) INTEGER */
/*          The number of base BETA digits in the mantissa of a */
/*          floating-point value. */

/*  EMIN    (input) INTEGER */
/*          The minimum exponent before (gradual) underflow. */

/*  IEEE    (input) LOGICAL */
/*          A logical flag specifying whether or not the arithmetic */
/*          system is thought to comply with the IEEE standard. */

/*  EMAX    (output) INTEGER */
/*          The largest exponent before overflow */

/*  RMAX    (output) DOUBLE PRECISION */
/*          The largest machine floating-point number. */

/* ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     First compute LEXP and UEXP, two powers of 2 that bound */
/*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum */
/*     approximately to the bound that is closest to abs(EMIN). */
/*     (EMAX is the exponent of the required number RMAX). */

    lexp = 1;
    exbits = 1;
L10:
    try_ = lexp << 1;
    if (try_ <= -(*emin)) {
	lexp = try_;
	++exbits;
	goto L10;
    }
    if (lexp == -(*emin)) {
	uexp = lexp;
    } else {
	uexp = try_;
	++exbits;
    }

/*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater */
/*     than or equal to EMIN. EXBITS is the number of bits needed to */
/*     store the exponent. */

    if (uexp + *emin > -lexp - *emin) {
	expsum = lexp << 1;
    } else {
	expsum = uexp << 1;
    }

/*     EXPSUM is the exponent range, approximately equal to */
/*     EMAX - EMIN + 1 . */

    *emax = expsum + *emin - 1;
    nbits = exbits + 1 + *p;

/*     NBITS is the total number of bits needed to store a */
/*     floating-point number. */

    if (nbits % 2 == 1 && *beta == 2) {

/*        Either there are an odd number of bits used to store a */
/*        floating-point number, which is unlikely, or some bits are 
*/
/*        not used in the representation of numbers, which is possible
, */
/*        (e.g. Cray machines) or the mantissa has an implicit bit, */
/*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
 */
/*        most likely. We have to assume the last alternative. */
/*        If this is true, then we need to reduce EMAX by one because 
*/
/*        there must be some way of representing zero in an implicit-b
it */
/*        system. On machines like Cray, we are reducing EMAX by one 
*/
/*        unnecessarily. */

	--(*emax);
    }

    if (*ieee) {

/*        Assume we are on an IEEE machine which reserves one exponent
 */
/*        for infinity and NaN. */

	--(*emax);
    }

/*     Now create RMAX, the largest machine number, which should */
/*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX . */

/*     First compute 1.0 - BETA**(-P), being careful that the */
/*     result is less than 1.0 . */

    recbas = 1. / *beta;
    z = *beta - 1.;
    y = 0.;
    i__1 = *p;
    for (i = 1; i <= i__1; ++i) {
	z *= recbas;
	if (y < 1.) {
	    oldy = y;
	}
	y = dlamc3_(&y, &z);
/* L20: */
    }
    if (y >= 1.) {
	y = oldy;
    }

/*     Now multiply by BETA**EMAX to get RMAX. */

    i__1 = *emax;
    for (i = 1; i <= i__1; ++i) {
	d__1 = y * *beta;
	y = dlamc3_(&d__1, &c_b21);
/* L30: */
    }

    *rmax = y;
    return 0;

/*     End of DLAMC5 */

} /* dlamc5_ */


logical lsame_(ca, cb, ca_len, cb_len)
char *ca, *cb;
ftnlen ca_len;
ftnlen cb_len;
{
    /* System generated locals */
    logical ret_val;

    /* Local variables */
    static integer inta, intb, zcode;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  LSAME returns .TRUE. if CA is the same letter as CB regardless of */
/*  case. */

/*  Arguments */
/*  ========= */

/*  CA      (input) CHARACTER*1 */
/*  CB      (input) CHARACTER*1 */
/*          CA and CB specify the single characters to be compared. */

/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test if the characters are equal */

    ret_val = *ca == *cb;
    if (ret_val) {
	return ret_val;
    }

/*     Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime */
/*     machines, on which ICHAR returns a value with bit 8 set. */
/*     ICHAR('A') on Prime machines returns 193 which is the same as */
/*     ICHAR('A') on an EBCDIC machine. */

    inta = *ca;
    intb = *cb;

    if (zcode == 90 || zcode == 122) {

/*        ASCII is assumed - ZCODE is the ASCII code of either lower o
r */
/*        upper case 'Z'. */

	if (inta >= 97 && inta <= 122) {
	    inta += -32;
	}
	if (intb >= 97 && intb <= 122) {
	    intb += -32;
	}

    } else if (zcode == 233 || zcode == 169) {

/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower
 or */
/*        upper case 'Z'. */

	if ((inta >= 129 && inta <= 137) || (inta >= 145 && inta <= 153) || (inta 
		>= 162 && inta <= 169)) {
	    inta += 64;
	}
	if ((intb >= 129 && intb <= 137) || (intb >= 145 && intb <= 153) || (intb 
		>= 162 && intb <= 169)) {
	    intb += 64;
	}

    } else if (zcode == 218 || zcode == 250) {

/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII cod
e */
/*        plus 128 of either lower or upper case 'Z'. */

	if (inta >= 225 && inta <= 250) {
	    inta += -32;
	}
	if (intb >= 225 && intb <= 250) {
	    intb += -32;
	}
    }
    ret_val = inta == intb;

/*     RETURN */

/*     End of LSAME */

    return ret_val;
} /* lsame_ */



/********************************* dbals.c ***********************************/
/* dblas.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Table of constant values */



/* Subroutine */ int daxpy_(n, da, dx, incx, dy, incy)
integer *n;
doublereal *da, *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     constant times a vector plus a vector. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dy[i] += *da * dx[i];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 4) {
	dy[i] += *da * dx[i];
	dy[i + 1] += *da * dx[i + 1];
	dy[i + 2] += *da * dx[i + 2];
	dy[i + 3] += *da * dx[i + 3];
/* L50: */
    }
    return 0;
} /* daxpy_ */

/* Subroutine */ int dcopy_(n, dx, incx, dy, incy)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     copies a vector, x, to a vector, y. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dy[iy] = dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dy[i] = dx[i];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 7) {
	dy[i] = dx[i];
	dy[i + 1] = dx[i + 1];
	dy[i + 2] = dx[i + 2];
	dy[i + 3] = dx[i + 3];
	dy[i + 4] = dx[i + 4];
	dy[i + 5] = dx[i + 5];
	dy[i + 6] = dx[i + 6];
/* L50: */
    }
    return 0;
} /* dcopy_ */

doublereal ddot_(n, dx, incx, dy, incy)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i, m;
    static doublereal dtemp;
    static integer ix, iy, mp1;


/*     forms the dot product of two vectors. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dtemp += dx[ix] * dy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dtemp += dx[i] * dy[i];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 5) {
	dtemp = dtemp + dx[i] * dy[i] + dx[i + 1] * dy[i + 1] + dx[i + 2] * 
		dy[i + 2] + dx[i + 3] * dy[i + 3] + dx[i + 4] * dy[i + 4];
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* ddot_ */






doublereal dnrm2_(n, dx, incx)
integer *n;
doublereal *dx;
integer *incx;
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal cutlo = 8.232e-11;
    static doublereal cuthi = 1.304e19;

    /* Format strings */
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_110[] = "";

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal xmax;
    static integer next, i, j, ix;
    static doublereal hitest, sum;

    /* Assigned format variables */
    char *next_fmt;

    /* Parameter adjustments */
    --dx;

    /* Function Body */

/*     euclidean norm of the n-vector stored in dx() with storage */
/*     increment incx . */
/*     if    n .le. 0 return with result = 0. */
/*     if n .ge. 1 then incx must be .ge. 1 */

/*           c.l.lawson, 1978 jan 08 */
/*     modified to correct failure to update ix, 1/25/92. */
/*     modified 3/93 to return if incx .le. 0. */

/*     four phase method     using two built-in constants that are */
/*     hopefully applicable to all machines. */
/*         cutlo = maximum of  dsqrt(u/eps)  over all known machines. */
/*         cuthi = minimum of  dsqrt(v)      over all known machines. */
/*     where */
/*         eps = smallest no. such that eps + 1. .gt. 1. */
/*         u   = smallest positive no.   (underflow limit) */
/*         v   = largest  no.            (overflow  limit) */

/*     brief outline of algorithm.. */

/*     phase 1    scans zero components. */
/*     move to phase 2 when a component is nonzero and .le. cutlo */
/*     move to phase 3 when a component is .gt. cutlo */
/*     move to phase 4 when a component is .ge. cuthi/m */
/*     where m = n for x() real and m = 2*n for complex. */

/*     values for cutlo and cuthi.. */
/*     from the environmental parameters listed in the imsl converter */
/*     document the limiting values are as follows.. */
/*     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are 
*/
/*                   univac and dec at 2**(-103) */
/*                   thus cutlo = 2**(-51) = 4.44089e-16 */
/*     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec. */
/*                   thus cuthi = 2**(63.5) = 1.30438e19 */
/*     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec. */
/*                   thus cutlo = 2**(-33.5) = 8.23181d-11 */
/*     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19 */
/*     data cutlo, cuthi / 8.232d-11,  1.304d19 / */
/*     data cutlo, cuthi / 4.441e-16,  1.304e19 / */

    if (*n > 0 && *incx > 0) {
	goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
    next_fmt = fmt_30;
    sum = zero;
    i = 1;
    ix = 1;
/*                                                 begin main loop */
L20:
    switch ((int)next) {
	case 0: goto L30;
	case 1: goto L50;
	case 2: goto L70;
	case 3: goto L110;
    }
L30:
    if ((d__1 = dx[i], abs(d__1)) > cutlo) {
	goto L85;
    }
    next = 1;
    next_fmt = fmt_50;
    xmax = zero;

/*                        phase 1.  sum is zero */

L50:
    if (dx[i] == zero) {
	goto L200;
    }
    if ((d__1 = dx[i], abs(d__1)) > cutlo) {
	goto L85;
    }

/*                                prepare for phase 2. */
    next = 2;
    next_fmt = fmt_70;
    goto L105;

/*                                prepare for phase 4. */

L100:
    ix = j;
    next = 3;
    next_fmt = fmt_110;
    sum = sum / dx[i] / dx[i];
L105:
    xmax = (d__1 = dx[i], abs(d__1));
    goto L115;

/*                   phase 2.  sum is small. */
/*                             scale to avoid destructive underflow. */

L70:
    if ((d__1 = dx[i], abs(d__1)) > cutlo) {
	goto L75;
    }

/*                     common code for phases 2 and 4. */
/*                     in phase 4 sum is large.  scale to avoid overflow. 
*/

L110:
    if ((d__1 = dx[i], abs(d__1)) <= xmax) {
	goto L115;
    }
/* Computing 2nd power */
    d__1 = xmax / dx[i];
    sum = one + sum * (d__1 * d__1);
    xmax = (d__1 = dx[i], abs(d__1));
    goto L200;

L115:
/* Computing 2nd power */
    d__1 = dx[i] / xmax;
    sum += d__1 * d__1;
    goto L200;


/*                  prepare for phase 3. */

L75:
    sum = sum * xmax * xmax;


/*     for real or d.p. set hitest = cuthi/n */
/*     for complex      set hitest = cuthi/(2*n) */

L85:
    hitest = cuthi / (real) (*n);

/*                   phase 3.  sum is mid-range.  no scaling. */

    i__1 = *n;
    for (j = ix; j <= i__1; ++j) {
	if ((d__1 = dx[i], abs(d__1)) >= hitest) {
	    goto L100;
	}
/* Computing 2nd power */
	d__1 = dx[i];
	sum += d__1 * d__1;
	i += *incx;
/* L95: */
    }
    ret_val = sqrt(sum);
    goto L300;

L200:
    ++ix;
    i += *incx;
    if (ix <= *n) {
	goto L20;
    }

/*              end of main loop. */

/*              compute square root and adjust for scaling. */

    ret_val = xmax * sqrt(sum);
L300:
    return ret_val;
} /* dnrm2_ */




/****************************** lib77.a **************************************/
/***************************** pow_di.c **************************************/

#ifdef KR_headers
double pow_di(ap, bp) doublereal *ap; integer *bp;
#else
double pow_di(doublereal *ap, integer *bp)
#endif
{
double pow, x;
integer n;

pow = 1;
x = *ap;
n = *bp;

if(n != 0)
	{
	if(n < 0)
		{
		n = -n;
		x = 1/x;
		}
	for( ; ; )
		{
		if(n & 01)
			pow *= x;
		if(n >>= 1)
			x *= x;
		else
			break;
		}
	}
return(pow);
}



/***************************** d_stop.c **************************************/
#undef abs




/***************************** sig_die.c *************************************/

#include "stdio.h"
#include "signal.h"

#ifndef SIGIOT
#define SIGIOT SIGABRT
#endif

#ifdef KR_headers
void sig_die(s, kill) register char *s; int kill;
#else
#include "stdlib.h"
#ifdef __cplusplus
extern "C" {
#endif
 extern void f_exit(void);

void sig_die(register char *s, int kill)
#endif
{
	/* print error message, then clear buffers */
	fprintf(stderr, "%s\n", s);
	fflush(stderr);
	f_exit();
	fflush(stderr);

	if(kill)
		{
		/* now get a core */
		signal(SIGIOT, SIG_DFL);
		abort();
		}
	else
		exit(1);
	}
#ifdef __cplusplus
}
#endif


/***************************** libI77.a **************************************/


#ifdef KR_headers
integer f_clos(a) cllist *a;
#else
#undef abs
#undef min
#undef max
#include "stdlib.h"
#ifdef NON_UNIX_STDIO
#ifndef unlink
#define unlink remove
#endif
#else
#ifdef MSDOS
#include "io.h"
#else
#ifdef __cplusplus
extern "C" int unlink(const char*);
#else
extern int unlink(const char*);
#endif
#endif
#endif

integer f_clos(cllist *a)
#endif
{	unit *b;

	if(a->cunit >= MXUNIT) return(0);
	b= &f__units[a->cunit];
	if(b->ufd==NULL)
		goto done;
	if (!a->csta) {
	  if (b->uscrtch == 1){ goto Delete;}
	  else		      { goto Keep; }
	}
	switch(*a->csta) {
		default:
	 	Keep:
		case 'k':
		case 'K':
			if(b->uwrt == 1)
				t_runc((alist *)a);
			if(b->ufnm) {
				fclose(b->ufd);
				free(b->ufnm);
				}
			break;
		case 'd':
		case 'D':
		Delete:
			if(b->ufnm) {
				fclose(b->ufd);
				unlink(b->ufnm); /*SYSDEP*/
				free(b->ufnm);
				}
		}
	b->ufd=NULL;
 done:
	b->uend=0;
	b->ufnm=NULL;
	return(0);
	}
 void
#ifdef KR_headers
f_exit()
#else
f_exit(void)
#endif
{	int i;
	static cllist xx;
	if (!xx.cerr) {
		xx.cerr=1;
		xx.csta=NULL;
		for(i=0;i<MXUNIT;i++)
		{
			xx.cunit=i;
			(void) f_clos(&xx);
		}
	}
}



/****************************** err.c ****************************************/

#ifndef NON_UNIX_STDIO
#include "sys/types.h"
#include "sys/stat.h"
#endif

#ifdef NON_UNIX_STDIO
#ifdef KR_headers
extern char *malloc();
#else
#undef abs
#undef min
#undef max
#include "stdlib.h"
#endif
#endif

/*global definitions*/
static unit f__units[MXUNIT];	/*unit table*/
static flag f__init;	/*0 on entry, 1 after initializations*/
static cilist *f__elist;	/*active external io list*/
static flag f__reading;	/*1 if reading, 0 if writing*/
static flag f__cplus,f__cblank;
static char *f__fmtbuf;
static flag f__external;	/*1 if external io, 0 if internal */
#ifdef KR_headers
int (*f__doed)(),(*f__doned)();
static int (*f__doend)(),(*f__donewrec)(),(*f__dorevert)();
static int (*f__getn)(),(*f__putn)();	/*for formatted io*/
#else
static int (*f__putn)(int);	/*for formatted io*/
static int (*f__doed)(struct f__syl*, char*, ftnlen),(*f__doned)(struct f__syl*);
static int (*f__dorevert)(void),(*f__donewrec)(void),(*f__doend)(void);
#endif
static flag f__sequential;	/*1 if sequential io, 0 if direct*/
static flag f__formatted;	/*1 if formatted io, 0 if unformatted*/
static FILE *f__cf;	/*current file*/
static unit *f__curunit;	/*current unit*/
static int f__recpos;	/*place in current record*/
static int f__cursor,f__scale;
#define MAXERR (sizeof(F_err)/sizeof(char *)+100)

int 
#ifdef KR_headers
f__canseek(f) FILE *f; /*SYSDEP*/
#else
f__canseek(FILE *f) /*SYSDEP*/
#endif
{
#ifdef NON_UNIX_STDIO
	return !isatty(fileno(f));
#else
	struct stat x;

	if (fstat(fileno(f),&x) < 0)
		return(0);
#ifdef S_IFMT
	switch(x.st_mode & S_IFMT) {
	case S_IFDIR:
	case S_IFREG:
		if(x.st_nlink > 0)	/* !pipe */
			return(1);
		else
			return(0);
	case S_IFCHR:
		if(isatty(fileno(f)))
			return(0);
		return(1);
#ifdef S_IFBLK
	case S_IFBLK:
		return(1);
#endif
	}
#else
#ifdef S_ISDIR
	/* POSIX version */
	if (S_ISREG(x.st_mode) || S_ISDIR(x.st_mode)) {
		if(x.st_nlink > 0)	/* !pipe */
			return(1);
		else
			return(0);
		}
	if (S_ISCHR(x.st_mode)) {
		if(isatty(fileno(f)))
			return(0);
		return(1);
		}
	if (S_ISBLK(x.st_mode))
		return(1);
#else
	Help! How does fstat work on this system?
#endif
#endif
	return(0);	/* who knows what it is? */
#endif
}

 void
#ifdef KR_headers
f__fatal(n,s) char *s;
#else
f__fatal(int n, char *s)
#endif
{
	if(n<100 && n>=0) perror(s); /*SYSDEP*/
	else if(n >= (int)MAXERR || n < -1)
	{	fprintf(stderr,"%s: illegal error number %d\n",s,n);
	}
	else if(n == -1) fprintf(stderr,"%s: end of file\n",s);
	else
		fprintf(stderr,"%s: %s\n",s,F_err[n-100]);
	if (f__curunit) {
		fprintf(stderr,"apparent state: unit %d ",f__curunit-f__units);
		fprintf(stderr, f__curunit->ufnm ? "named %s\n" : "(unnamed)\n",
			f__curunit->ufnm);
		}
	else
		fprintf(stderr,"apparent state: internal I/O\n");
	if (f__fmtbuf)
		fprintf(stderr,"last format: %s\n",f__fmtbuf);
	fprintf(stderr,"lately %s %s %s %s",f__reading?"reading":"writing",
		f__sequential?"sequential":"direct",f__formatted?"formatted":"unformatted",
		f__external?"external":"internal");
	sig_die(" IO", 1);
}
/*initialization routine*/
 VOID
f_init(Void)
{	unit *p;

	f__init=1;
	p= &f__units[0];
	p->ufd=stderr;
	p->useek=f__canseek(stderr);
#ifdef NON_UNIX_STDIO
	setbuf(stderr, (char *)malloc(BUFSIZ));
#else
	//stderr->_flag &= ~_IONBF;
#endif
	p->ufmt=1;
	p->uwrt=1;
	p = &f__units[5];
	p->ufd=stdin;
	p->useek=f__canseek(stdin);
	p->ufmt=1;
	p->uwrt=0;
	p= &f__units[6];
	p->ufd=stdout;
	p->useek=f__canseek(stdout);
	/* IOLBUF and setvbuf only in system 5+ */
#ifdef COMMENTED_OUT
	if(isatty(fileno(stdout))) {
		extern char _sobuf[];
		setbuf(stdout, _sobuf);
		/* setvbuf(stdout, _IOLBF, 0, 0); */
		   /* the buf arg in setvbuf? */
		/* p->useek = 1; */
		/* only within a record no bigger than BUFSIZ */
	}
#endif
	p->ufmt=1;
	p->uwrt=1;
}


int 
#ifdef KR_headers
f__nowwriting(x) unit *x;
#else
f__nowwriting(unit *x)
#endif
{
	long loc;
	int k;


	if (!x->ufnm)
		goto cantwrite;
	if (x->uwrt == 3) { /* just did write, rewind */
#ifdef NON_UNIX_STDIO
		if (!(f__cf = x->ufd =
				freopen(x->ufnm,f__w_mode[x->ufmt],x->ufd)))
#else
		if (close(creat(x->ufnm,0666)))
#endif
			goto cantwrite;
		}
	else {
		loc=ftell(x->ufd);
#ifdef NON_UNIX_STDIO
		if (!(f__cf = x->ufd =
			freopen(x->ufnm, f__w_mode[x->ufmt+2], x->ufd)))
#else
		if (fclose(x->ufd) < 0
		|| (k = x->uwrt == 2 ? creat(x->ufnm,0666)
				     : open(x->ufnm,O_WRONLY)) < 0
		|| (f__cf = x->ufd = fdopen(k,f__w_mode[x->ufmt])) == NULL)
#endif
			{
			x->ufd = NULL;
 cantwrite:
			errno = 127;
			return(1);
			}
		(void) fseek(x->ufd,loc,SEEK_SET);
		}
	x->uwrt = 1;
	return(0);
}

 int
#ifdef KR_headers
err__fl(f, m, s) int f, m; char *s;
#else
err__fl(int f, int m, char *s)
#endif
{
	if (!f)
		f__fatal(m, s);
	if (f__doend)
		(*f__doend)();
	return errno = m;
	}


/****************************** err.c ****************************************/


#include "sys/types.h"


#ifdef KR_headers
extern char *strcpy();
#else
#undef abs
#undef min
#undef max
#include "stdlib.h"
#include "string.h"
#endif

#ifdef NON_UNIX_STDIO
#ifndef unlink
#define unlink remove
#endif
#else
#ifdef MSDOS
#include "io.h"
#endif
#endif

#ifdef NON_UNIX_STDIO
extern char *f__r_mode[], *f__w_mode[];
#endif


 static int
#ifdef NON_UNIX_STDIO
#ifdef KR_headers
copy(from, len, to) char *from, *to; register long len;
#else
copy(FILE *from, register long len, FILE *to)
#endif
{
	int k, len1;
	char buf[BUFSIZ];

	while(fread(buf, len1 = len > BUFSIZ ? BUFSIZ : (int)len, 1, from)) {
		if (!fwrite(buf, len1, 1, to))
			return 1;
		if ((len -= len1) <= 0)
			break;
		}
	return 0;
	}
#else
#ifdef KR_headers
copy(from, len, to) char *from, *to; register long len;
#else
copy(char *from, register long len, char *to)
#endif
{
	register int n;
	int k, rc = 0, tmp;
	char buf[BUFSIZ];

	if ((k = open(from, O_RDONLY)) < 0)
		return 1;
	if ((tmp = creat(to,0666)) < 0)
		return 1;
	while((n = read(k, buf, len > BUFSIZ ? BUFSIZ : (int)len)) > 0) {
		if (write(tmp, buf, n) != n)
			{ rc = 1; break; }
		if ((len -= n) <= 0)
			break;
		}
	close(k);
	close(tmp);
	return n < 0 ? 1 : rc;
	}
#endif

#ifndef L_tmpnam
#define L_tmpnam 16
#endif

 int
#ifdef KR_headers
t_runc(a) alist *a;
#else
t_runc(alist *a)
#endif
{
	char nm[L_tmpnam+12];	/* extra space in case L_tmpnam is tiny */
	long loc, len;
	unit *b;
#ifdef NON_UNIX_STDIO
	FILE *bf, *tf;
#else
	FILE *bf;
#endif
	int rc = 0;

	b = &f__units[a->aunit];
	if(b->url)
		return(0);	/*don't truncate direct files*/
	loc=ftell(bf = b->ufd);
	fseek(bf,0L,SEEK_END);
	len=ftell(bf);
	if (loc >= len || b->useek == 0 || b->ufnm == NULL)
		return(0);
#ifdef NON_UNIX_STDIO
	fclose(b->ufd);
#else
	rewind(b->ufd);	/* empty buffer */
#endif
	if (!loc) {
#ifdef NON_UNIX_STDIO
		if (!(bf = fopen(b->ufnm, f__w_mode[b->ufmt])))
#else
		if (close(creat(b->ufnm,0666)))
#endif
			rc = 1;
		if (b->uwrt)
			b->uwrt = 1;
		goto done;
		}
#ifdef _POSIX_SOURCE
	tmpnam(nm);
#else
		{
		  static int n = 100000;
		  sprintf(nm,"tmp.F%d",n++);
		  if ( n >= 499999 ) n = 100000;
		  //strcpy(nm,"tmp.FXXXXXX");
		  //mktemp(nm);
		}

#endif
#ifdef NON_UNIX_STDIO
	if (!(bf = fopen(b->ufnm, f__r_mode[0]))) {
 bad:
		rc = 1;
		goto done;
		}
	if (!(tf = fopen(nm, f__w_mode[0])))
		goto bad;
	if (copy(bf, loc, tf)) {
 bad1:
		rc = 1;
		goto done1;
		}
	if (!(bf = freopen(b->ufnm, f__w_mode[0], bf)))
		goto bad1;
	if (!(tf = freopen(nm, f__r_mode[0], tf)))
		goto bad1;
	if (copy(tf, loc, bf))
		goto bad1;
	if (f__w_mode[0] != f__w_mode[b->ufmt]) {
	 	if (!(bf = freopen(b->ufnm, f__w_mode[b->ufmt+2], bf)))
			goto bad1;
		fseek(bf, loc, SEEK_SET);
		}
done1:
	fclose(tf);
	unlink(nm);
done:
	f__cf = b->ufd = bf;
#else
	if (copy(b->ufnm, loc, nm)
	 || copy(nm, loc, b->ufnm))
		rc = 1;
	unlink(nm);
	fseek(b->ufd, loc, SEEK_SET);
done:
#endif
	if (rc)
		err(a->aerr,111,"endfile");
	return 0;
	}


/****************************** open.c ***************************************/
#include "sys/types.h"
#include "sys/stat.h"
#include "string.h"

#ifdef KR_headers
extern char *malloc(), *mktemp();
extern integer f_clos();
#else
#undef abs
#undef min
#undef max
#include "stdlib.h"
extern int f__canseek(FILE*);
extern integer f_clos(cllist*);
#endif

int 
#ifdef KR_headers
f__isdev(s) char *s;
#else
f__isdev(char *s)
#endif
{
#ifdef NON_UNIX_STDIO
	int i, j;

	i = open(s,O_RDONLY);
	if (i == -1)
		return 0;
	j = isatty(i);
	close(i);
	return j;
#else
	struct stat x;

	if(stat(s, &x) == -1) return(0);
#ifdef S_IFMT
	switch(x.st_mode&S_IFMT) {
		case S_IFREG:
		case S_IFDIR:
			return(0);
		}
#else
#ifdef S_ISREG
	/* POSIX version */
	if(S_ISREG(x.st_mode) || S_ISDIR(x.st_mode))
		return(0);
	else
#else
	Help! How does stat work on this system?
#endif
#endif
		return(1);
#endif
}
#ifdef KR_headers
integer f_open(a) olist *a;
#else
integer f_open(olist *a)
#endif
{	unit *b;
	int n;
	integer rv;
	char buf[256], *s;
	cllist x;
#ifdef NON_UNIX_STDIO
	FILE *tf;
#else
	struct stat stb;
#endif
	if(a->ounit>=MXUNIT || a->ounit<0)
		err(a->oerr,101,"open")
	f__curunit = b = &f__units[a->ounit];
	if(b->ufd) {
		if(a->ofnm==0)
		{
		same:	if (a->oblnk)
				b->ublnk = *a->oblnk == 'z' || *a->oblnk == 'Z';
			return(0);
		}
#ifdef NON_UNIX_STDIO
		if (b->ufnm
		 && strlen(b->ufnm) == a->ofnmlen
		 && !strncmp(b->ufnm, b->ufnm, (unsigned)a->ofnmlen))
			goto same;
#else
		g_char(a->ofnm,a->ofnmlen,buf);
		if (f__inode(buf,&n) == b->uinode && n == b->udev)
			goto same;
#endif
		x.cunit=a->ounit;
		x.csta=0;
		x.cerr=a->oerr;
		if ((rv = f_clos(&x)) != 0)
			return rv;
		}
	b->url = (int)a->orl;
	b->ublnk = a->oblnk && (*a->oblnk == 'z' || *a->oblnk == 'Z');
	if(a->ofm==0)
	{	if(b->url>0) b->ufmt=0;
		else b->ufmt=1;
	}
	else if(*a->ofm=='f' || *a->ofm == 'F') b->ufmt=1;
	else b->ufmt=0;
#ifdef url_Adjust
	if (b->url && !b->ufmt)
		url_Adjust(b->url);
#endif
	if (a->ofnm) {
		g_char(a->ofnm,a->ofnmlen,buf);
		if (!buf[0])
			err(a->oerr,107,"open")
		}
	else
		sprintf(buf, "fort.%ld", a->ounit);
	b->uscrtch = 0;
	switch(a->osta ? *a->osta : 'u')
	{
	case 'o':
	case 'O':
#ifdef NON_UNIX_STDIO
		if(access(buf,0))
#else
		if(stat(buf,&stb))
#endif
			err(a->oerr,errno,"open")
		break;
	 case 's':
	 case 'S':
		b->uscrtch=1;
#ifdef _POSIX_SOURCE
		tmpnam(buf);
#else
		{
		  static int n = 500000;
		  sprintf(buf,"tmp.F%d",n++);
		  if ( n >= 999999 ) n = 500000;
		  //(void) strcpy(buf,"tmp.FXXXXXX");
		  //(void) mktemp(buf);
		}
#endif
		goto replace;
	case 'n':
	case 'N':
#ifdef NON_UNIX_STDIO
		if(!access(buf,0))
#else
		if(!stat(buf,&stb))
#endif
			err(a->oerr,128,"open")
		/* no break */
	case 'r':	/* Fortran 90 replace option */
	case 'R':
 replace:
#ifdef NON_UNIX_STDIO
		if (tf = fopen(buf,f__w_mode[0]))
			fclose(tf);
#else
		(void) close(creat(buf, 0666));
#endif
	}

	b->ufnm=(char *) malloc((unsigned int)(strlen(buf)+1));
	if(b->ufnm==NULL) err(a->oerr,113,"no space");
	(void) strcpy(b->ufnm,buf);
	b->uend=0;
	b->uwrt = 0;
	if(f__isdev(buf))
	{	b->ufd = fopen(buf,f__r_mode[b->ufmt]);
		if(b->ufd==NULL) err(a->oerr,errno,buf)
	}
	else {
		if(!(b->ufd = fopen(buf, f__r_mode[b->ufmt]))) {
#ifdef NON_UNIX_STDIO
			if (b->ufd = fopen(buf, f__w_mode[b->ufmt+2]))
				b->uwrt = 2;
			else if (b->ufd = fopen(buf, f__w_mode[b->ufmt]))
				b->uwrt = 1;
			else
#else
			if ((n = open(buf,O_WRONLY)) >= 0)
				b->uwrt = 2;
			else {
				n = creat(buf, 0666);
				b->uwrt = 1;
				}
			if (n < 0
			|| (b->ufd = fdopen(n, f__w_mode[b->ufmt])) == NULL)
#endif
				err(a->oerr, errno, "open");
			}
	}
	b->useek=f__canseek(b->ufd);
#ifndef NON_UNIX_STDIO
	if((b->uinode=f__inode(buf,&b->udev))==-1)
		err(a->oerr,108,"open")
#endif
       if(b->useek) {
	  if (a->orl) {
			rewind(b->ufd);
	  }
	  else {if ((s = a->oacc) && (*s == 'a' || *s == 'A')
			&& fseek(b->ufd, 0L, SEEK_END))
				err(a->oerr,129,"open");
	  }
       }
	return(0);
}


int 
#ifdef KR_headers
fk_open(seq,fmt,n) ftnint n;
#else
fk_open(int seq, int fmt, ftnint n)
#endif
{	char nbuf[10];
	olist a;
	(void) sprintf(nbuf,"fort.%ld",n);
	a.oerr=1;
	a.ounit=n;
	a.ofnm=nbuf;
	a.ofnmlen=strlen(nbuf);
	a.osta=NULL;
	a.oacc= seq==SEQ?"s":"d";
	a.ofm = fmt==FMT?"f":"u";
	a.orl = seq==DIR?1:0;
	a.oblnk=NULL;
	return(f_open(&a));
}

/****************************** util.c ***************************************/
#ifndef NON_UNIX_STDIO
#include "sys/types.h"
#include "sys/stat.h"
#endif

 VOID
#ifdef KR_headers
g_char(a,alen,b) char *a,*b; ftnlen alen;
#else
g_char(char *a, ftnlen alen, char *b)
#endif
{
	char *x = a + alen, *y = b + alen;

	for(;; y--) {
		if (x <= a) {
			*b = 0;
			return;
			}
		if (*--x != ' ')
			break;
		}
	*y-- = 0;
	do *y-- = *x;
		while(x-- > a);
}

#ifndef NON_UNIX_STDIO
#ifdef KR_headers
long f__inode(a, dev) char *a; int *dev;
#else
long f__inode(char *a, int *dev)
#endif
{	struct stat x;
	if(stat(a,&x)<0) return(-1);
	*dev = x.st_dev;
	return(x.st_ino);
}
#endif


/****************************** wsfe.c ***************************************/

/*write sequential formatted external*/
extern int f__hiwater;

int 
#ifdef KR_headers
x_putc(c)
#else
x_putc(int c)
#endif
{
	/* this uses \n as an indicator of record-end */
	if(c == '\n' && f__recpos < f__hiwater) {	/* fseek calls fflush, a loss */
#ifndef NON_UNIX_STDIO
	  /*
		if(f__cf->_ptr + f__hiwater - f__recpos < buf_end(f__cf))
			f__cf->_ptr += f__hiwater - f__recpos;
		else
	  */
#endif
			(void) fseek(f__cf, (long)(f__hiwater - f__recpos), SEEK_CUR);
	}
	f__recpos++;
	return putc(c,f__cf);
}
int x_wSL(Void)
{
	(*f__putn)('\n');
	f__recpos=0;
	f__cursor = 0;
	f__hiwater = 0;
	return(1);
}
int xw_end(Void)
{
	if(f__nonl == 0)
		(*f__putn)('\n');
	f__hiwater = f__recpos = f__cursor = 0;
	return(0);
}
int xw_rev(Void)
{
	if(f__workdone) (*f__putn)('\n');
	f__hiwater = f__recpos = f__cursor = 0;
	return(f__workdone=0);
}

#ifdef KR_headers
integer s_wsfe(a) cilist *a;	/*start*/
#else
integer s_wsfe(cilist *a)	/*start*/
#endif
{	int n;
	if(!f__init) f_init();
	if((n=c_sfe(a))&&(n!=0)) return(n);
	f__reading=0;
	f__sequential=1;
	f__formatted=1;
	f__external=1;
	f__elist=a;
	f__hiwater = f__cursor=f__recpos=0;
	f__nonl = 0;
	f__scale=0;
	f__fmtbuf=a->cifmt;
	f__curunit = &f__units[a->ciunit];
	f__cf=f__curunit->ufd;
	if(pars_f(f__fmtbuf)<0) err(a->cierr,100,"startio");
	f__putn= x_putc;
	f__doed= w_ed;
	f__doned= w_ned;
	f__doend=xw_end;
	f__dorevert=xw_rev;
	f__donewrec=x_wSL;
	fmt_bg();
	f__cplus=0;
	f__cblank=f__curunit->ublnk;
	if(f__curunit->uwrt != 1 && f__nowwriting(f__curunit))
		err(a->cierr,errno,"write start");
	return(0);
}

/****************************** wrtfmt.c *************************************/

extern int f__cursor;
#ifdef KR_headers
extern char *f__icvt();
#else
extern char *f__icvt(long, int*, int*, int);
#endif
static int f__hiwater;
static icilist *f__svic;
static char *f__icptr;
int 
mv_cur(Void)	/* shouldn't use fseek because it insists on calling fflush */
		/* instead we know too much about stdio */
{
	if(f__external == 0) {
		if(f__cursor < 0) {
			if(f__hiwater < f__recpos)
				f__hiwater = f__recpos;
			f__recpos += f__cursor;
			f__icptr += f__cursor;
			f__cursor = 0;
			if(f__recpos < 0)
				err(f__elist->cierr, 110, "left off");
		}
		else if(f__cursor > 0) {
			if(f__recpos + f__cursor >= f__svic->icirlen)
				err(f__elist->cierr, 110, "recend");
			if(f__hiwater <= f__recpos)
				for(; f__cursor > 0; f__cursor--)
					(*f__putn)(' ');
			else if(f__hiwater <= f__recpos + f__cursor) {
				f__cursor -= f__hiwater - f__recpos;
				f__icptr += f__hiwater - f__recpos;
				f__recpos = f__hiwater;
				for(; f__cursor > 0; f__cursor--)
					(*f__putn)(' ');
			}
			else {
				f__icptr += f__cursor;
				f__recpos += f__cursor;
			}
			f__cursor = 0;
		}
		return(0);
	}
	if(f__cursor > 0) {
		if(f__hiwater <= f__recpos)
			for(;f__cursor>0;f__cursor--) (*f__putn)(' ');
		else if(f__hiwater <= f__recpos + f__cursor) {
#ifndef NON_UNIX_STDIO
		  /*
		  if(f__cf->_ptr + f__hiwater - f__recpos < buf_end(f__cf))
				f__cf->_ptr += f__hiwater - f__recpos;
			else
		  */
#endif
				(void) fseek(f__cf, (long) (f__hiwater - f__recpos), SEEK_CUR);
			f__cursor -= f__hiwater - f__recpos;
			f__recpos = f__hiwater;
			for(; f__cursor > 0; f__cursor--)
				(*f__putn)(' ');
		}
		else {
#ifndef NON_UNIX_STDIO
		  /*
			if(f__cf->_ptr + f__cursor < buf_end(f__cf))
				f__cf->_ptr += f__cursor;
			else
		  */
#endif
				(void) fseek(f__cf, (long)f__cursor, SEEK_CUR);
			f__recpos += f__cursor;
		}
	}
	if(f__cursor<0)
	{
		if(f__cursor+f__recpos<0) err(f__elist->cierr,110,"left off");
#ifndef NON_UNIX_STDIO
		/*
		if(f__cf->_ptr + f__cursor >= f__cf->_base)
			f__cf->_ptr += f__cursor;
		else
		*/
#endif
		if(f__curunit && f__curunit->useek)
			(void) fseek(f__cf,(long)f__cursor,SEEK_CUR);
		else
			err(f__elist->cierr,106,"fmt");
		if(f__hiwater < f__recpos)
			f__hiwater = f__recpos;
		f__recpos += f__cursor;
		f__cursor=0;
	}
	return(0);
}

 static int
#ifdef KR_headers
wrt_Z(n,w,minlen,len) Uint *n; int w, minlen; ftnlen len;
#else
wrt_Z(Uint *n, int w, int minlen, ftnlen len)
#endif
{
	register char *s, *se;

	register int i, w1;
	static int one = 1;
	static char hex[] = "0123456789ABCDEF";
	s = (char *)n;
	--len;
	if (*(char *)&one) {
		/* little endian */
		se = s;
		s += len;
		i = -1;
		}
	else {
		se = s + len;
		i = 1;
		}
	for(;; s += i)
		if (s == se || *s)
			break;
	w1 = (i*(se-s) << 1) + 1;
	if (*s & 0xf0)
		w1++;
	if (w1 > w)
		for(i = 0; i < w; i++)
			(*f__putn)('*');
	else {
		if ((minlen -= w1) > 0)
			w1 += minlen;
		while(--w >= w1)
			(*f__putn)(' ');
		while(--minlen >= 0)
			(*f__putn)('0');
		if (!(*s & 0xf0)) {
			(*f__putn)(hex[*s & 0xf]);
			if (s == se)
				return 0;
			s += i;
			}
		for(;; s += i) {
			(*f__putn)(hex[*s >> 4 & 0xf]);
			(*f__putn)(hex[*s & 0xf]);
			if (s == se)
				break;
			}
		}
	return 0;
	}

 static int
#ifdef KR_headers
wrt_I(n,w,len, base) Uint *n; ftnlen len; register int base;
#else
wrt_I(Uint *n, int w, ftnlen len, register int base)
#endif
{	int ndigit,sign,spare,i;
	long x;
	char *ans;
	if(len==sizeof(integer)) x=n->il;
	else if(len == sizeof(char)) x = n->ic;
#ifdef Allow_TYQUAD
	else if (len == sizeof(longint)) x = n->ili;
#endif
	else x=n->is;
	ans=f__icvt(x,&ndigit,&sign, base);
	spare=w-ndigit;
	if(sign || f__cplus) spare--;
	if(spare<0)
		for(i=0;i<w;i++) (*f__putn)('*');
	else
	{	for(i=0;i<spare;i++) (*f__putn)(' ');
		if(sign) (*f__putn)('-');
		else if(f__cplus) (*f__putn)('+');
		for(i=0;i<ndigit;i++) (*f__putn)(*ans++);
	}
	return(0);
}
 static int
#ifdef KR_headers
wrt_IM(n,w,m,len,base) Uint *n; ftnlen len; int base;
#else
wrt_IM(Uint *n, int w, int m, ftnlen len, int base)
#endif
{	int ndigit,sign,spare,i,xsign;
	long x;
	char *ans;
	if(sizeof(integer)==len) x=n->il;
	else if(len == sizeof(char)) x = n->ic;
	else x=n->is;
	ans=f__icvt(x,&ndigit,&sign, base);
	if(sign || f__cplus) xsign=1;
	else xsign=0;
	if(ndigit+xsign>w || m+xsign>w)
	{	for(i=0;i<w;i++) (*f__putn)('*');
		return(0);
	}
	if(x==0 && m==0)
	{	for(i=0;i<w;i++) (*f__putn)(' ');
		return(0);
	}
	if(ndigit>=m)
		spare=w-ndigit-xsign;
	else
		spare=w-m-xsign;
	for(i=0;i<spare;i++) (*f__putn)(' ');
	if(sign) (*f__putn)('-');
	else if(f__cplus) (*f__putn)('+');
	for(i=0;i<m-ndigit;i++) (*f__putn)('0');
	for(i=0;i<ndigit;i++) (*f__putn)(*ans++);
	return(0);
}
 static int
#ifdef KR_headers
wrt_AP(s) char *s;
#else
wrt_AP(char *s)
#endif
{	char quote;
	if(f__cursor && mv_cur()) return(mv_cur());
	quote = *s++;
	for(;*s;s++)
	{	if(*s!=quote) (*f__putn)(*s);
		else if(*++s==quote) (*f__putn)(*s);
		else return(1);
	}
	return(1);
}
 static int
#ifdef KR_headers
wrt_H(a,s) char *s;
#else
wrt_H(int a, char *s)
#endif
{
	if(f__cursor && mv_cur()) return(mv_cur());
	while(a--) (*f__putn)(*s++);
	return(1);
}
int 
#ifdef KR_headers
wrt_L(n,len, sz) Uint *n; ftnlen sz;
#else
wrt_L(Uint *n, int len, ftnlen sz)
#endif
{	int i;
	long x;
	if(sizeof(long)==sz) x=n->il;
	else if(sz == sizeof(char)) x = n->ic;
	else x=n->is;
	for(i=0;i<len-1;i++)
		(*f__putn)(' ');
	if(x) (*f__putn)('T');
	else (*f__putn)('F');
	return(0);
}
 static int
#ifdef KR_headers
wrt_A(p,len) char *p; ftnlen len;
#else
wrt_A(char *p, ftnlen len)
#endif
{
	while(len-- > 0) (*f__putn)(*p++);
	return(0);
}
 static int
#ifdef KR_headers
wrt_AW(p,w,len) char * p; ftnlen len;
#else
wrt_AW(char * p, int w, ftnlen len)
#endif
{
	while(w>len)
	{	w--;
		(*f__putn)(' ');
	}
	while(w-- > 0)
		(*f__putn)(*p++);
	return(0);
}

 static int
#ifdef KR_headers
wrt_G(p,w,d,e,len) ufloat *p; ftnlen len;
#else
wrt_G(ufloat *p, int w, int d, int e, ftnlen len)
#endif
{	double up = 1,x;
	int i,oldscale=f__scale,n,j;
	x= len==sizeof(real)?p->pf:p->pd;
	if(x < 0 ) x = -x;
	if(x<.1) return(wrt_E(p,w,d,e,len));
	for(i=0;i<=d;i++,up*=10)
	{	if(x>=up) continue;
		f__scale=0;
		if(e==0) n=4;
		else	n=e+2;
		i=wrt_F(p,w-n,d-i,len);
		for(j=0;j<n;j++) (*f__putn)(' ');
		f__scale=oldscale;
		return(i);
	}
	return(wrt_E(p,w,d,e,len));
}

int 
#ifdef KR_headers
w_ed(p,ptr,len) struct f__syl *p; char *ptr; ftnlen len;
#else
w_ed(struct f__syl *p, char *ptr, ftnlen len)
#endif
{
	if(f__cursor && mv_cur()) return(mv_cur());
	switch(p->op)
	{
	default:
		fprintf(stderr,"w_ed, unexpected code: %d\n", p->op);
		sig_die(f__fmtbuf, 1);
	case I:	return(wrt_I((Uint *)ptr,p->p1,len, 10));
	case IM:
		return(wrt_IM((Uint *)ptr,p->p1,p->p2,len,10));

		/* O and OM don't work right for character, double, complex, */
		/* or doublecomplex, and they differ from Fortran 90 in */
		/* showing a minus sign for negative values. */

	case O:	return(wrt_I((Uint *)ptr, p->p1, len, 8));
	case OM:
		return(wrt_IM((Uint *)ptr,p->p1,p->p2,len,8));
	case L:	return(wrt_L((Uint *)ptr,p->p1, len));
	case A: return(wrt_A(ptr,len));
	case AW:
		return(wrt_AW(ptr,p->p1,len));
	case D:
	case E:
	case EE:
		return(wrt_E((ufloat *)ptr,p->p1,p->p2,p->p3,len));
	case G:
	case GE:
		return(wrt_G((ufloat *)ptr,p->p1,p->p2,p->p3,len));
	case F:	return(wrt_F((ufloat *)ptr,p->p1,p->p2,len));

		/* Z and ZM assume 8-bit bytes. */

	case Z: return(wrt_Z((Uint *)ptr,p->p1,0,len));
	case ZM:
		return(wrt_Z((Uint *)ptr,p->p1,p->p2,len));
	}
}
int 
#ifdef KR_headers
w_ned(p) struct f__syl *p;
#else
w_ned(struct f__syl *p)
#endif
{
	switch(p->op)
	{
	default: fprintf(stderr,"w_ned, unexpected code: %d\n", p->op);
		sig_die(f__fmtbuf, 1);
	case SLASH:
		return((*f__donewrec)());
	case T: f__cursor = p->p1-f__recpos - 1;
		return(1);
	case TL: f__cursor -= p->p1;
		if(f__cursor < -f__recpos)	/* TL1000, 1X */
			f__cursor = -f__recpos;
		return(1);
	case TR:
	case X:
		f__cursor += p->p1;
		return(1);
	case APOS:
		return(wrt_AP(*(char **)&p->p2));
	case H:
		return(wrt_H(p->p1,*(char **)&p->p2));
	}
}

/****************************** fmt.c ****************************************/
#define skip(s) while(*s==' ') s++
#ifdef interdata
#define SYLMX 300
#endif
#ifdef pdp11
#define SYLMX 300
#endif
#ifdef vax
#define SYLMX 300
#endif
#ifndef SYLMX
#define SYLMX 300
#endif
#define GLITCH '\2'
	/* special quote character for stu */
extern int f__cursor,f__scale;
extern flag f__cblank,f__cplus;	/*blanks in I and compulsory plus*/
static struct f__syl f__syl[SYLMX];
static int f__parenlvl,f__pc,f__revloc;

#ifdef KR_headers
static char *ap_end(s) char *s;
#else
static char *ap_end(char *s)
#endif
{	char quote;
	quote= *s++;
	for(;*s;s++)
	{	if(*s!=quote) continue;
		if(*++s!=quote) return(s);
	}
	if(f__elist->cierr) {
		errno = 100;
		return(NULL);
	}
	f__fatal(100, "bad string");
	/*NOTREACHED*/ return 0;
}
int 
#ifdef KR_headers
op_gen(a,b,c,d)
#else
op_gen(int a, int b, int c, int d)
#endif
{	struct f__syl *p= &f__syl[f__pc];
	if(f__pc>=SYLMX)
	{	fprintf(stderr,"format too complicated:\n");
		sig_die(f__fmtbuf, 1);
	}
	p->op=a;
	p->p1=b;
	p->p2=c;
	p->p3=d;
	return(f__pc++);
}
#ifdef KR_headers
char *f_list();
char *gt_num(s,n) char *s; int *n;
#else
char *f_list(char*);
char *gt_num(char *s, int *n)
#endif
{	int m=0,f__cnt=0;
	char c;
	for(c= *s;;c = *s)
	{	if(c==' ')
		{	s++;
			continue;
		}
		if(c>'9' || c<'0') break;
		m=10*m+c-'0';
		f__cnt++;
		s++;
	}
	if(f__cnt==0) *n=1;
	else *n=m;
	return(s);
}
#ifdef KR_headers
char *f_s(s,curloc) char *s;
#else
char *f_s(char *s, int curloc)
#endif
{
	skip(s);
	if(*s++!='(')
	{
		return(NULL);
	}
	if(f__parenlvl++ ==1) f__revloc=curloc;
	if(op_gen(RET1,curloc,0,0)<0 ||
		(s=f_list(s))==NULL)
	{
		return(NULL);
	}
	skip(s);
	return(s);
}
int 
#ifdef KR_headers
ne_d(s,p) char *s,**p;
#else
ne_d(char *s, char **p)
#endif
{	int n,x,sign=0;
	struct f__syl *sp;
	switch(*s)
	{
	default:
		return(0);
	case ':': (void) op_gen(COLON,0,0,0); break;
	case '$':
		(void) op_gen(NONL, 0, 0, 0); break;
	case 'B':
	case 'b':
		if(*++s=='z' || *s == 'Z') (void) op_gen(BZ,0,0,0);
		else (void) op_gen(BN,0,0,0);
		break;
	case 'S':
	case 's':
		if(*(s+1)=='s' || *(s+1) == 'S')
		{	x=SS;
			s++;
		}
		else if(*(s+1)=='p' || *(s+1) == 'P')
		{	x=SP;
			s++;
		}
		else x=S;
		(void) op_gen(x,0,0,0);
		break;
	case '/': (void) op_gen(SLASH,0,0,0); break;
	case '-': sign=1;
	case '+':	s++;	/*OUTRAGEOUS CODING TRICK*/
	case '0': case '1': case '2': case '3': case '4':
	case '5': case '6': case '7': case '8': case '9':
		s=gt_num(s,&n);
		switch(*s)
		{
		default:
			return(0);
		case 'P':
		case 'p': if(sign) n= -n; (void) op_gen(P,n,0,0); break;
		case 'X':
		case 'x': (void) op_gen(X,n,0,0); break;
		case 'H':
		case 'h':
			sp = &f__syl[op_gen(H,n,0,0)];
			*(char **)&sp->p2 = s + 1;
			s+=n;
			break;
		}
		break;
	case GLITCH:
	case '"':
	case '\'':
		sp = &f__syl[op_gen(APOS,0,0,0)];
		*(char **)&sp->p2 = s;
		if((*p = ap_end(s)) == NULL)
			return(0);
		return(1);
	case 'T':
	case 't':
		if(*(s+1)=='l' || *(s+1) == 'L')
		{	x=TL;
			s++;
		}
		else if(*(s+1)=='r'|| *(s+1) == 'R')
		{	x=TR;
			s++;
		}
		else x=T;
		s=gt_num(s+1,&n);
		s--;
		(void) op_gen(x,n,0,0);
		break;
	case 'X':
	case 'x': (void) op_gen(X,1,0,0); break;
	case 'P':
	case 'p': (void) op_gen(P,1,0,0); break;
	}
	s++;
	*p=s;
	return(1);
}
int 
#ifdef KR_headers
e_d(s,p) char *s,**p;
#else
e_d(char *s, char **p)
#endif
{	int i,im,n,w,d,e,found=0,x=0;
	char *sv=s;
	s=gt_num(s,&n);
	(void) op_gen(STACK,n,0,0);
	switch(*s++)
	{
	default: break;
	case 'E':
	case 'e':	x=1;
	case 'G':
	case 'g':
		found=1;
		s=gt_num(s,&w);
		if(w==0) break;
		if(*s=='.')
		{	s++;
			s=gt_num(s,&d);
		}
		else d=0;
		if(*s!='E' && *s != 'e')
			(void) op_gen(x==1?E:G,w,d,0);	/* default is Ew.dE2 */
		else
		{	s++;
			s=gt_num(s,&e);
			(void) op_gen(x==1?EE:GE,w,d,e);
		}
		break;
	case 'O':
	case 'o':
		i = O;
		im = OM;
		goto finish_I;
	case 'Z':
	case 'z':
		i = Z;
		im = ZM;
		goto finish_I;
	case 'L':
	case 'l':
		found=1;
		s=gt_num(s,&w);
		if(w==0) break;
		(void) op_gen(L,w,0,0);
		break;
	case 'A':
	case 'a':
		found=1;
		skip(s);
		if(*s>='0' && *s<='9')
		{	s=gt_num(s,&w);
			if(w==0) break;
			(void) op_gen(AW,w,0,0);
			break;
		}
		(void) op_gen(A,0,0,0);
		break;
	case 'F':
	case 'f':
		found=1;
		s=gt_num(s,&w);
		if(w==0) break;
		if(*s=='.')
		{	s++;
			s=gt_num(s,&d);
		}
		else d=0;
		(void) op_gen(F,w,d,0);
		break;
	case 'D':
	case 'd':
		found=1;
		s=gt_num(s,&w);
		if(w==0) break;
		if(*s=='.')
		{	s++;
			s=gt_num(s,&d);
		}
		else d=0;
		(void) op_gen(D,w,d,0);
		break;
	case 'I':
	case 'i':
		i = I;
		im = IM;
 finish_I:
		found=1;
		s=gt_num(s,&w);
		if(w==0) break;
		if(*s!='.')
		{	(void) op_gen(i,w,0,0);
			break;
		}
		s++;
		s=gt_num(s,&d);
		(void) op_gen(im,w,d,0);
		break;
	}
	if(found==0)
	{	f__pc--; /*unSTACK*/
		*p=sv;
		return(0);
	}
	*p=s;
	return(1);
}
#ifdef KR_headers
char *i_tem(s) char *s;
#else
char *i_tem(char *s)
#endif
{	char *t;
	int n,curloc;
	if(*s==')') return(s);
	if(ne_d(s,&t)) return(t);
	if(e_d(s,&t)) return(t);
	s=gt_num(s,&n);
	if((curloc=op_gen(STACK,n,0,0))<0) return(NULL);
	return(f_s(s,curloc));
}
#ifdef KR_headers
char *f_list(s) char *s;
#else
char *f_list(char *s)
#endif
{
	for(;*s!=0;)
	{	skip(s);
		if((s=i_tem(s))==NULL) return(NULL);
		skip(s);
		if(*s==',') s++;
		else if(*s==')')
		{	if(--f__parenlvl==0)
			{
				(void) op_gen(REVERT,f__revloc,0,0);
				return(++s);
			}
			(void) op_gen(GOTO,0,0,0);
			return(++s);
		}
	}
	return(NULL);
}

int 
#ifdef KR_headers
pars_f(s) char *s;
#else
pars_f(char *s)
#endif
{
	f__parenlvl=f__revloc=f__pc=0;
	if(f_s(s,0) == NULL)
	{
		return(-1);
	}
	return(0);
}
#define STKSZ 10
static int f__cnt[STKSZ],f__ret[STKSZ],f__cp,f__rp;
static flag f__workdone, f__nonl;

int 
#ifdef KR_headers
type_f(n)
#else
type_f(int n)
#endif
{
	switch(n)
	{
	default:
		return(n);
	case RET1:
		return(RET1);
	case REVERT: return(REVERT);
	case GOTO: return(GOTO);
	case STACK: return(STACK);
	case X:
	case SLASH:
	case APOS: case H:
	case T: case TL: case TR:
		return(NED);
	case F:
	case I:
	case IM:
	case A: case AW:
	case O: case OM:
	case L:
	case E: case EE: case D:
	case G: case GE:
	case Z: case ZM:
		return(ED);
	}
}
#ifdef KR_headers
integer do_fio(number,ptr,len) ftnint *number; ftnlen len; char *ptr;
#else
integer do_fio(ftnint *number, char *ptr, ftnlen len)
#endif
{	struct f__syl *p;
	int n,i;
	for(i=0;i<*number;i++,ptr+=len)
	{
loop:	switch(type_f((p= &f__syl[f__pc])->op))
	{
	default:
		fprintf(stderr,"unknown code in do_fio: %d\n%s\n",
			p->op,f__fmtbuf);
		err(f__elist->cierr,100,"do_fio");
	case NED:
		if((*f__doned)(p))
		{	f__pc++;
			goto loop;
		}
		f__pc++;
		continue;
	case ED:
		if(f__cnt[f__cp]<=0)
		{	f__cp--;
			f__pc++;
			goto loop;
		}
		if(ptr==NULL)
			return((*f__doend)());
		f__cnt[f__cp]--;
		f__workdone=1;
		if((n=(*f__doed)(p,ptr,len))>0)
			errfl(f__elist->cierr,errno,"fmt");
		if(n<0)
			err(f__elist->ciend,(EOF),"fmt");
		continue;
	case STACK:
		f__cnt[++f__cp]=p->p1;
		f__pc++;
		goto loop;
	case RET1:
		f__ret[++f__rp]=p->p1;
		f__pc++;
		goto loop;
	case GOTO:
		if(--f__cnt[f__cp]<=0)
		{	f__cp--;
			f__rp--;
			f__pc++;
			goto loop;
		}
		f__pc=1+f__ret[f__rp--];
		goto loop;
	case REVERT:
		f__rp=f__cp=0;
		f__pc = p->p1;
		if(ptr==NULL)
			return((*f__doend)());
		if(!f__workdone) return(0);
		if((n=(*f__dorevert)()) != 0) return(n);
		goto loop;
	case COLON:
		if(ptr==NULL)
			return((*f__doend)());
		f__pc++;
		goto loop;
	case NONL:
		f__nonl = 1;
		f__pc++;
		goto loop;
	case S:
	case SS:
		f__cplus=0;
		f__pc++;
		goto loop;
	case SP:
		f__cplus = 1;
		f__pc++;
		goto loop;
	case P:	f__scale=p->p1;
		f__pc++;
		goto loop;
	case BN:
		f__cblank=0;
		f__pc++;
		goto loop;
	case BZ:
		f__cblank=1;
		f__pc++;
		goto loop;
	}
	}
	return(0);
}
int 
en_fio(Void)
{	ftnint one=1;
	return(do_fio(&one,(char *)NULL,(ftnint)0));
}
 VOID
fmt_bg(Void)
{
	f__workdone=f__cp=f__rp=f__pc=f__cursor=0;
	f__cnt[0]=f__ret[0]=0;
}


/****************************** sfe.c ****************************************/
/* sequential formatted external common routines*/

extern char *f__fmtbuf;

integer e_rsfe(Void)
{	int n;
	n=en_fio();
	if (f__cf == stdout)
		fflush(stdout);
	else if (f__cf == stderr)
		fflush(stderr);
	f__fmtbuf=NULL;
	return(n);
}
#ifdef KR_headers
static int c_sfe(a) cilist *a; /* check */
#else
static int c_sfe(cilist *a) /* check */
#endif
{	unit *p;
	if(a->ciunit >= MXUNIT || a->ciunit<0)
		err(a->cierr,101,"startio");
	p = &f__units[a->ciunit];
	if(p->ufd==NULL && fk_open(SEQ,FMT,a->ciunit)) err(a->cierr,114,"sfe")
	if(!p->ufmt) err(a->cierr,102,"sfe")
	return(0);
}
integer e_wsfe(Void)
{	return(e_rsfe());
}

/****************************** fmtlib.c *************************************/
/*	@(#)fmtlib.c	1.2	*/
#define MAXINTLENGTH 23
#ifdef KR_headers
char *f__icvt(value,ndigit,sign, base) long value; int *ndigit,*sign;
 register int base;
#else
char *f__icvt(long value, int *ndigit, int *sign, int base)
#endif
{	static char buf[MAXINTLENGTH+1];
	register int i;
	if(value>0) *sign=0;
	else if(value<0)
	{	value = -value;
		*sign= 1;
	}
	else
	{	*sign=0;
		*ndigit=1;
		buf[MAXINTLENGTH]='0';
		return(&buf[MAXINTLENGTH]);
	}
	for(i=MAXINTLENGTH-1;value>0;i--)
	{	*(buf+i)=(int)(value%base)+'0';
		value /= base;
	}
	*ndigit=MAXINTLENGTH-1-i;
	return(&buf[i+1]);
}


/****************************** wref.c ***************************************/

#ifndef VAX
#include "ctype.h"
#endif

#ifndef KR_headers
#undef abs
#undef min
#undef max
#include "stdlib.h"
#include "string.h"
#endif

int 
#ifdef KR_headers
wrt_E(p,w,d,e,len) ufloat *p; ftnlen len;
#else
wrt_E(ufloat *p, int w, int d, int e, ftnlen len)
#endif
{
	char buf[FMAX+EXPMAXDIGS+4], *s, *se;
	int d1, delta, e1, i, sign, signspace;
	double dd;
#ifndef VAX
	int e0 = e;
#endif

	if(e <= 0)
		e = 2;
	if(f__scale) {
		if(f__scale >= d + 2 || f__scale <= -d)
			goto nogood;
		}
	if(f__scale <= 0)
		--d;
	if (len == sizeof(real))
		dd = p->pf;
	else
		dd = p->pd;
	if (dd < 0.) {
		signspace = sign = 1;
		dd = -dd;
		}
	else {
		sign = 0;
		signspace = (int)f__cplus;
#ifndef VAX
		if (!dd)
			dd = 0.;	/* avoid -0 */
#endif
		}
	delta = w - (2 /* for the . and the d adjustment above */
			+ 2 /* for the E+ */ + signspace + d + e);
	if (delta < 0) {
nogood:
		while(--w >= 0)
			PUT('*');
		return(0);
		}
	if (f__scale < 0)
		d += f__scale;
	if (d > FMAX) {
		d1 = d - FMAX;
		d = FMAX;
		}
	else
		d1 = 0;
	sprintf(buf,"%#.*E", d, dd);
#ifndef VAX
	/* check for NaN, Infinity */
	if (!isdigit(buf[0])) {
		switch(buf[0]) {
			case 'n':
			case 'N':
				signspace = 0;	/* no sign for NaNs */
			}
		delta = w - strlen(buf) - signspace;
		if (delta < 0)
			goto nogood;
		while(--delta >= 0)
			PUT(' ');
		if (signspace)
			PUT(sign ? '-' : '+');
		for(s = buf; *s; s++)
			PUT(*s);
		return 0;
		}
#endif
	se = buf + d + 3;
	if (f__scale != 1 && dd)
		sprintf(se, "%+.2d", atoi(se) + 1 - f__scale);
	s = ++se;
	if (e < 2) {
		if (*s != '0')
			goto nogood;
		}
#ifndef VAX
	/* accommodate 3 significant digits in exponent */
	if (s[2]) {
#ifdef Pedantic
		if (!e0 && !s[3])
			for(s -= 2, e1 = 2; s[0] = s[1]; s++);

	/* Pedantic gives the behavior that Fortran 77 specifies,	*/
	/* i.e., requires that E be specified for exponent fields	*/
	/* of more than 3 digits.  With Pedantic undefined, we get	*/
	/* the behavior that Cray displays -- you get a bigger		*/
	/* exponent field if it fits.	*/
#else
		if (!e0) {
			for(s -= 2, e1 = 2; (s[0] = s[1]) && (s[0] !=0); s++)
#ifdef CRAY
				delta--;
			if ((delta += 4) < 0)
				goto nogood
#endif
				;
			}
#endif
		else if (e0 >= 0)
			goto shift;
		else
			e1 = e;
		}
	else
 shift:
#endif
		for(s += 2, e1 = 2; *s; ++e1, ++s)
			if (e1 >= e)
				goto nogood;
	while(--delta >= 0)
		PUT(' ');
	if (signspace)
		PUT(sign ? '-' : '+');
	s = buf;
	i = f__scale;
	if (f__scale <= 0) {
		PUT('.');
		for(; i < 0; ++i)
			PUT('0');
		PUT(*s);
		s += 2;
		}
	else if (f__scale > 1) {
		PUT(*s);
		s += 2;
		while(--i > 0)
			PUT(*s++);
		PUT('.');
		}
	if (d1) {
		se -= 2;
		while(s < se) PUT(*s++);
		se += 2;
		do PUT('0'); while(--d1 > 0);
		}
	while(s < se)
		PUT(*s++);
	if (e < 2)
		PUT(s[1]);
	else {
		while(++e1 <= e)
			PUT('0');
		while(*s)
			PUT(*s++);
		}
	return 0;
	}

int 
#ifdef KR_headers
wrt_F(p,w,d,len) ufloat *p; ftnlen len;
#else
wrt_F(ufloat *p, int w, int d, ftnlen len)
#endif
{
	int d1, sign, n;
	double x;
	char *b, buf[MAXINTDIGS+MAXFRACDIGS+4], *s;

	x= (len==sizeof(real)?p->pf:p->pd);
	if (d < MAXFRACDIGS)
		d1 = 0;
	else {
		d1 = d - MAXFRACDIGS;
		d = MAXFRACDIGS;
		}
	if (x < 0.)
		{ x = -x; sign = 1; }
	else {
		sign = 0;
#ifndef VAX
		if (!x)
			x = 0.;
#endif
		}

	 n = f__scale; {
		if (n > 0)
			do x *= 10.; while(--n > 0);
		else
			do x *= 0.1; while(++n < 0);
	}
		
#ifdef USE_STRLEN
	sprintf(b = buf, "%#.*f", d, x);
	n = strlen(b) + d1;
#else
	n = sprintf(b = buf, "%#.*f", d, x) + d1;
#endif

	if (buf[0] == '0' && d)
		{ ++b; --n; }
	if (sign) {
		/* check for all zeros */
		for(s = b;;) {
			while(*s == '0') s++;
			switch(*s) {
				case '.':
					s++; continue;
				case 0:
					sign = 0;
				}
			break;
			}
		}
	if (sign || f__cplus)
		++n;
	if (n > w) {
		while(--w >= 0)
			PUT('*');
		return 0;
		}
	for(w -= n; --w >= 0; )
		PUT(' ');
	if (sign)
		PUT('-');
	else if (f__cplus)
		PUT('+');
	while((n = *b++) && (n != 0))
		PUT(n);
	while(--d1 >= 0)
		PUT('0');
	return 0;
	}

