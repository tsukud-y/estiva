/* 
 * PSC98 input data generation
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define MATMAX	4000000
#define NZMAX	9

/* input data generatoin */

static	int	key = -1;

static int yosen1();
static int yosen2();
static int yosen3();
static int yosen4();
static int yosen5();
static int sample();

static void genmat_common(i, ja, a, b)
    int i, *ja;
    double *a, *b;
{
	char	*buff;

	if ( key == -1 ) {
		buff = getenv( "PSC98" );
		if ( buff == NULL ) buff = "0";
		key = atoi(buff);
	}

	switch (key) {
		case 1: yosen1(i, ja, a, b); break;
		case 2: yosen2(i, ja, a, b); break;
		case 3: yosen3(i, ja, a, b); break;
		case 4: yosen4(i, ja, a, b); break;
		case 5: yosen5(i, ja, a, b); break;
		default: key=0; sample(i, ja, a, b); break;
	}
}

/*
 *  Fortran77 Interface Routine
 */
/*
 *  call genmat(i, ja, a, b)
 */
void
genmat_(i, ja, a, b)
    int    *i, ja[];
    double a[], *b;
{
#if IDX_BEGINS_WITH_ZERO
    genmat_common(*i, ja, a, b);
#else
    int i_prime;
    if (*i == -1) i_prime = -1;
    else i_prime = *i - 1;
    genmat_common(i_prime, ja, a, b);
#endif
}

/*
 *  C Interface Routine
 */

void genmat(i, ja, a, b)
    int i, *ja;
    double *a, *b;
{
#if IDX_BEGINS_WITH_ZERO
    genmat_common(i, ja, a, b);
#else
    int i_prime;
    if (i == -1) i_prime = -1;
    else i_prime = i - 1;
    genmat_common(i_prime, ja, a, b);
#endif
}


void
chkval(s, n, x)
    FILE * s;			/* output stream */
    int n;
    double *x;
{
    double a[NZMAX];
    int ja[NZMAX], i, mat_size;
    double b, resmax, max_r = 0.0, l_time;

    genmat(-1, ja, a, &b);
    mat_size = ja[0];
    resmax = b;

    assert(n == mat_size);

    for (i = 0; i < mat_size; ++i) {
	double r;
	int j = 0;

#if IDX_BEGINS_WITH_ZERO
	genmat(i, ja, a, &b);
#else
	genmat(i + 1, ja, a, &b);
#endif
	r = b;

	for (j = 0; ja[j] != -1 && j < NZMAX; ++j) {
#if IDX_BEGINS_WITH_ZERO
	    r -= a[j] * x[ja[j]];
#else
	    r -= a[j] * x[ja[j]-1];
#endif
	}
	r = fabs(r);

	if (!(max_r >= r)) {
	  if (max_r >= 0) {
	    max_r = r;
	  }
	}
    }

    fprintf(s, "Problem NO : %d\n", key);
    fprintf(s, "|b - Ax|_inf = %g", max_r);

    if (max_r < resmax) {
	fprintf (s, "  (OK)\n");
    }
    else {
	fprintf (s, "  (NG)\n");
    }
}
/*
 *  Simple sort
 */

static void
sort_PSC98(n, idx, perm)
    int n, *idx, *perm;
{
    int i, j;

    for (i = 0; i < n; ++i)
	perm[i] = i;

    for (i = 0; i < n; ++i) {
	int min = idx[perm[i]];
	int min_idx = i;
	for (j = i + 1; j < n; ++j) {
	    if (min > idx[perm[j]]) {
		min = idx[perm[j]];
		min_idx = j;
	    }
	}

	{
	    int t = perm[i];
	    perm[i] = perm[min_idx];
	    perm[min_idx] = t;
	}
    }
}

#define SIZE_X		127
#define SIZE_Y		127
#define PROB_SIZE	(SIZE_X * SIZE_Y)
#define MAX_NZERO	5
#define RESMAX		1e-5

static int
sample(i, ja, a, b)
    int i, *ja;
    double *a, *b;
{
    double inv_h = (SIZE_X + 1) * (SIZE_X + 1);
    double d0 = 4.0 * inv_h;
    double d1 = -1.0 * inv_h;

    if (i == -1) {
	ja[0] = PROB_SIZE;
	ja[1] = PROB_SIZE
	    + 2 * (PROB_SIZE - SIZE_Y)
	    + 2 * (PROB_SIZE - SIZE_X);
	ja[2] = MAX_NZERO;
	*b = RESMAX;

	return 0;
    }
    else if (i >= 0 && i < PROB_SIZE) {
	int k = 0;

	if (i - SIZE_X >= 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - SIZE_X;
#else
	    ja[k] = i + 1 - SIZE_X;
#endif
	    a[k] = d1;
	    ++k;
	}
	if (i % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - 1;
#else
	    ja[k] = i;
#endif
	    a[k] = d1;
	    ++k;
	}
#if IDX_BEGINS_WITH_ZERO
	ja[k] = i;
#else
	ja[k] = i + 1;
#endif
	a[k] = d0;
	++k;
	if ((i + 1) % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + 1;
#else
	    ja[k] = i + 2;
#endif
	    a[k] = d1;
	    ++k;
	}
	if (i + SIZE_X < PROB_SIZE) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + SIZE_X;
#else
	    ja[k] = i + 1 + SIZE_X;
#endif
	    a[k] = d1;
	    ++k;
	}
	for (; k < NZMAX; ++k) {
	    ja[k] = -1;
	    a[k] = 0.0;
	}

	*b = 0.5 * (sin((double)i) + 1.0);

	return 0;
    }
    else {
	perror("sample");

	return 1;
    }
}

#undef SIZE_X
#undef SIZE_Y
#undef PROB_SIZE
#undef MAX_NZERO
#undef RESMAX

/************** yosen1 *************/

#define SIZE_X		511
#define SIZE_Y		511
#define PROB_SIZE	(SIZE_X * SIZE_Y)
#define MAX_NZERO	5
#define RESMAX		1e-5

static int
yosen1(i, ja, a, b)
    int i, *ja;
    double *a, *b;
{
    double inv_h = (SIZE_X + 1) * (SIZE_X + 1);
    double d0 = 4.0 * inv_h;
    double d1 = -1.0 * inv_h;

    if (i == -1) {
	ja[0] = PROB_SIZE;
	ja[1] = PROB_SIZE
	    + 2 * (PROB_SIZE - SIZE_Y)
	    + 2 * (PROB_SIZE - SIZE_X);
	ja[2] = MAX_NZERO;
	*b = RESMAX;

	return 0;
    }
    else if (i >= 0 && i < PROB_SIZE) {
	int k = 0;

	if (i - SIZE_X >= 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - SIZE_X;
#else
	    ja[k] = i + 1 - SIZE_X;
#endif
	    a[k] = d1;
	    ++k;
	}
	if (i % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - 1;
#else
	    ja[k] = i;
#endif
	    a[k] = d1;
	    ++k;
	}
#if IDX_BEGINS_WITH_ZERO
	ja[k] = i;
#else
	ja[k] = i + 1;
#endif
	a[k] = d0;
	++k;
	if ((i + 1) % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + 1;
#else
	    ja[k] = i + 2;
#endif
	    a[k] = d1;
	    ++k;
	}
	if (i + SIZE_X < PROB_SIZE) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + SIZE_X;
#else
	    ja[k] = i + 1 + SIZE_X;
#endif
	    a[k] = d1;
	    ++k;
	}
	for (; k < NZMAX; ++k) {
	    ja[k] = -1;
	    a[k] = 0.0;
	}

	*b = 0.5 * (sin((double)i) + 1.0);

	return 0;
    }
    else {
	perror("yosen1");

	return 1;
    }
}

#undef SIZE_X
#undef SIZE_Y
#undef PROB_SIZE
#undef MAX_NZERO
#undef RESMAX

/************** yosen2 *************/

/*
 *  2-dimensional Poisson Equation with Dirichlet Condition
 *
 *  ___________
 *  | SIZE_X2 |
 *  |         |
 *  | SIZE_X1   |  SIZE_Y
 *  |           |
 *   ~~~~~~~~~~~
 */


#define SIZE_X1		511
#define SIZE_X2		245
#define SIZE_Y		511
#define MAX_NZERO	5
#define	RESMAX		1e-5

static int
yosen2(i, ja, a, b)
    int i, *ja;
    double *a, *b;
{
    double inv_h = (SIZE_Y + 1) * (SIZE_Y + 1);
    double d0 = 4.0 * inv_h;
    double d1 = -1.0 * inv_h;
    int n2 = SIZE_Y / 2;
    int N2 = SIZE_X1 * n2;
    int N = N2 + SIZE_X2 * (SIZE_Y - n2);

    if (i == -1) {
	ja[0] = N;
	ja[1] = N + 2 * (N - SIZE_Y)
	    + N - SIZE_X1
	    + N - SIZE_X2
	    - abs(SIZE_X1 - SIZE_X2);
	ja[2] = MAX_NZERO;
	*b = RESMAX;

	return 0;
    }
    else if (i >= 0 && i < N) {
	int k = 0;

	if (i >= SIZE_X1
	    && (i >= N2 + SIZE_X2 || i < N2 + SIZE_X1)) {
	    if (i >= N2 + SIZE_X2) {
#if IDX_BEGINS_WITH_ZERO
		ja[k] = i - SIZE_X2;
#else
		ja[k] = i + 1 - SIZE_X2;
#endif
		a[k] = d1;
		++k;
	    }
	    else {
#if IDX_BEGINS_WITH_ZERO
		ja[k] = i - SIZE_X1;
#else
		ja[k] = i + 1 - SIZE_X1;
#endif
		a[k] = d1;
		++k;
	    }
	}

	if (i >= N2 && (i - N2) % SIZE_X2 > 0
	    || i < N2 && i % SIZE_X1 > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - 1;
#else
	    ja[k] = i;
#endif
	    a[k] = d1;
	    ++k;
	}

	/**  Diagonal elements  **/

#if IDX_BEGINS_WITH_ZERO
	ja[k] = i;
#else
	ja[k] = i + 1;
#endif
	a[k] = d0;
	++k;

	if (i >= N2 && (i + 1 - N2) % SIZE_X2 > 0
	    || i < N2 && (i + 1) % SIZE_X1 > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + 1;
#else
	    ja[k] = i + 2;
#endif
	    a[k] = d1;
	    ++k;
	}

	if (i < N - SIZE_X2
	    && (i >= N2 || i < N2 - SIZE_X1 + SIZE_X2)) {
	    if (i >= N2) {
#if IDX_BEGINS_WITH_ZERO
		ja[k] = i + SIZE_X2;
#else
		ja[k] = i + 1 + SIZE_X2;
#endif
		a[k] = d1;
		++k;
	    }
	    else {
#if IDX_BEGINS_WITH_ZERO
		ja[k] = i + SIZE_X1;
#else
		ja[k] = i + 1 + SIZE_X1;
#endif
		a[k] = d1;
		++k;
	    }
	}

	for (; k < NZMAX; ++k) {
	    ja[k] = -1;
	    a[k] = 0.0;
	}

	/*  right-hand side  */

	*b = 0.5 * (sin((double)i) + 1.0);

	return 0;
    }
    else {
	perror("yosen2");

	return 1;
    }
}
#undef SIZE_X1
#undef SIZE_X2
#undef SIZE_Y
#undef PROB_SIZE
#undef MAX_NZERO
#undef RESMAX

/************** yosen3 *************/

/*
 *  3-dimensional Poisson Equation with Dirichlet Condition
 *
 *  
 */

#define SIZE_X		63
#define SIZE_Y		127
#define SIZE_Z		127
#define PROB_SIZE	(SIZE_X * SIZE_Y * SIZE_Z)
#define SIZE_XY		(SIZE_X * SIZE_Y)
#define MAX_NZERO	7
#define	RESMAX		1e-5

static int
yosen3(i, ja, a, b)
    int i, *ja;
    double *a, *b;
{
    double inv_h = (SIZE_X + 1) * (SIZE_X + 1);
    double d0 = 6.0 * inv_h;
    double d1 = -1.0 * inv_h;

    if (i == -1) {
	ja[0] = PROB_SIZE;
	ja[1] = PROB_SIZE
	    + 2 * (PROB_SIZE - SIZE_XY)
	    + 2 * (PROB_SIZE - SIZE_Y * SIZE_Z)
	    + 2 * (PROB_SIZE - SIZE_Z * SIZE_X);
	ja[2] = MAX_NZERO;
	*b = RESMAX;

	return 0;
    }
    else if (i >= 0 && i < PROB_SIZE) {
	int k = 0;

	if (i >= SIZE_XY) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - SIZE_XY;
#else
	    ja[k] = i + 1 - SIZE_XY;
#endif
	    a[k] = d1;
	    ++k;
	}
	if (i % SIZE_XY >= SIZE_X) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - SIZE_X;
#else
	    ja[k] = i + 1 - SIZE_X;
#endif
	    a[k] = d1;
	    ++k;
	}
	if (i % SIZE_XY % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - 1;
#else
	    ja[k] = i;
#endif
	    a[k] = d1;
	    ++k;
	}
#if IDX_BEGINS_WITH_ZERO
	ja[k] = i;
#else
	ja[k] = i + 1;
#endif
	a[k] = d0;
	++k;
	if ((i + 1) % SIZE_XY % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + 1;
#else
	    ja[k] = i + 2;
#endif
	    a[k] = d1;
	    ++k;
	}
	if (i % SIZE_XY < SIZE_XY - SIZE_X) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + SIZE_X;
#else
	    ja[k] = i + 1 + SIZE_X;
#endif
	    a[k] = d1;
	    ++k;
	}
	if (i < PROB_SIZE - SIZE_XY) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + SIZE_XY;
#else
	    ja[k] = i + 1 + SIZE_XY;
#endif
	    a[k] = d1;
	    ++k;
	}
	for (; k < NZMAX; ++k) {
	    ja[k] = -1;
	    a[k] = 0.0;
	}

	*b = 0.5 * (sin((double)i) + 1.0);

	return 0;
    }
    else {
	perror("yosen3");

	return 1;
    }
}
#undef SIZE_X
#undef SIZE_Y
#undef SIZE_Z
#undef SIZE_XY
#undef PROB_SIZE
#undef MAX_NZERO
#undef RESMAX

/************** yosen4 *************/

/*
 *  2-dimensional Poisson Equation with Dirichlet Condition
 *
 *  Non-uniform diffusion constant
 */

#define SIZE_X		511
#define SIZE_Y		511
#define PROB_SIZE	(SIZE_X * SIZE_Y)
#define MAX_NZERO	5
#define RESMAX          1e-5

/*
 *  +----------+
 *  | k1       |
 *  |    +--+  |
 *  |    |k2|  |
 *  |    +--+  |
 *  |          |
 *  |          |
 *  +----------+
 */

static void
get_coeff(i, c)
    int i;
    double *c;
{
    double k1 = 1;
    double k2 = 100;
    int ix, iy;
    int cx1, cy1, cx2, cy2;

    iy = i / SIZE_X;
    ix = i % SIZE_X;
    cx1 = SIZE_X / 4;
    cy1 = SIZE_Y / 4;
    cx2 = 3 * SIZE_X / 4;
    cy2 = 3 * SIZE_Y / 4;

    if ((ix > cx1 && iy > cy1)
	&& (ix <= cx2 && iy <= cy2))
	c[0] = k2;
    else
	c[0] = k1;
	
    if ((ix >= cx1 && iy > cy1)
	&& (ix < cx2 && iy <= cy2))
	c[1] = k2;
    else
	c[1] = k1;

    if ((ix > cx1 && iy >= cy1)
	&& (ix <= cx2 && iy < cy2))
	c[2] = k2;
    else
	c[2] = k1;

    if ((ix >= cx1 && iy >= cy1)
	&& (ix < cx2 && iy < cy2))
	c[3] = k2;
    else
	c[3] = k1;
}

static int
yosen4(i, ja, a, b)
    int i, *ja;
    double *a, *b;
{
    double inv_h = (SIZE_X + 1) * (SIZE_X + 1);

    if (i == -1) {
	ja[0] = PROB_SIZE;
	ja[1] = PROB_SIZE
	    + 2 * (PROB_SIZE - SIZE_Y)
	    + 2 * (PROB_SIZE - SIZE_X);
	ja[2] = MAX_NZERO;
	*b = RESMAX;

	return 0;
    }
    else if (i >= 0 && i < PROB_SIZE) {
	int k = 0;
	double c[4];
	double c0, c1, c2, c3;

	/*
	 *  c[2] | c[3]
	 *  -----+-----
	 *  c[0] | c[1]
	 */

	get_coeff(i, c);

	c0 = c[0] * inv_h;
	c1 = c[1] * inv_h;
	c2 = c[2] * inv_h;
	c3 = c[3] * inv_h;

	if (i - SIZE_X >= 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - SIZE_X;
#else
	    ja[k] = i + 1 - SIZE_X;
#endif
	    a[k] = -.5 * (c0 + c1);
	    ++k;
	}
	if (i % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - 1;
#else
	    ja[k] = i;
#endif
	    a[k] = -.5 * (c0 + c2);
	    ++k;
	}
#if IDX_BEGINS_WITH_ZERO
	ja[k] = i;
#else
	ja[k] = i + 1;
#endif
	a[k] = c0 + c1 + c2 + c3;
	++k;
	if ((i + 1) % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + 1;
#else
	    ja[k] = i + 2;
#endif
	    a[k] = -.5 * (c1 + c3);
	    ++k;
	}
	if (i + SIZE_X < PROB_SIZE) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + SIZE_X;
#else
	    ja[k] = i + 1 + SIZE_X;
#endif
	    a[k] = -.5 * (c2 + c3);
	    ++k;
	}
	for (; k < NZMAX; ++k) {
	    ja[k] = -1;
	    a[k] = 0.0;
	}

	*b = 0.5 * (sin((double)i) + 1.0);

	return 0;
    }
    else {
	perror("genmat");

	return 1;
    }
}

#undef SIZE_X
#undef SIZE_Y
#undef PROB_SIZE
#undef MAX_NZERO
#undef RESMAX

/************** yosen5 *************/

/*
 *  1-dimensional Poisson Equation with Dirichlet Condition
 *
 *  
 */

/* #define SIZE_X		524287 */
#define SIZE_X	        16383
#define PROB_SIZE	(SIZE_X)
#define MAX_NZERO	3
#define	RESMAX		1e-5

static int
yosen5(i, ja, a, b)
    int i, *ja;
    double *a, *b;
{
    double inv_h = (double)(SIZE_X + 1) * (double)(SIZE_X + 1);
    double d0 = 2.0 * inv_h;
    double d1 = -1.0 * inv_h;

    if (i == -1) {
	ja[0] = PROB_SIZE;
	ja[1] = PROB_SIZE + 2 * (PROB_SIZE - 1);
	ja[2] = MAX_NZERO;
	*b = RESMAX;

	return 0;
    }
    else if (i >= 0 && i < PROB_SIZE) {
	int k = 0;

	if (i - 1 >= 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - 1;
#else
	    ja[k] = i;
#endif
	    a[k] = d1;
	    ++k;
	}
#if IDX_BEGINS_WITH_ZERO
	ja[k] = i;
#else
	ja[k] = i + 1;
#endif
	a[k] = d0;
	++k;
	if (i + 1 < PROB_SIZE) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + 1;
#else
	    ja[k] = i + 2;
#endif
	    a[k] = d1;
	    ++k;
	}
	for (; k < NZMAX; ++k) {
	    ja[k] = -1;
	    a[k] = 0.0;
	}

	*b = 0.5 * (sin((double)i) + 1.0);

	return 0;
    }
    else {
	perror("yosen5");

	return 1;
    }
}
#undef SIZE_X
#undef PROB_SIZE
#undef MAX_NZERO
#undef RESMAX

