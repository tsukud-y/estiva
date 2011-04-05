/* 
 * PSC98 HONSEN input data generation
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define MATMAX	2000000
#define NZMAX	9

/* input data generatoin */

static	int	key = -1;

static int honsen1();
static int honsen2();
static int honsen3();
static int honsen4();
static int honsen5();
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
		case 1: honsen1(i, ja, a, b); break;
		case 2: honsen2(i, ja, a, b); break;
		case 3: honsen3(i, ja, a, b); break;
		case 4: honsen4(i, ja, a, b); break;
		case 5: honsen5(i, ja, a, b); break;
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
    double b, resmax, max_r = 0.0;

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

/*
 *  Right-hand term
 */

static double
trand1(i)
    int i;
{
    unsigned int t = 2 * i + 1;
    t *= 48828125;
    t *= 48828125;
    t *= 48828125;
    return .5 * (double)t / (double)(unsigned int)0x80000000;
}

static double
trand2(i)
    int i;
{
    unsigned int t = 2 * i + 1;
    t *= 1812433253;
    t *= 1812433253;
    t *= 1812433253;
    return .5 * (double)t / (double)(unsigned int)0x80000000;
}

static double
trand3(i)
    int i;
{
    unsigned int t = 2 * i + 1;
    t *= 1566083941;
    t *= 1566083941;
    t *= 1566083941;
    return .5 * (double)t / (double)(unsigned int)0x80000000;
}

static double
trand4(i)
    int i;
{
    unsigned int t = 2 * i + 1;
    t *= 69069;
    t *= 69069;
    t *= 69069;
    t *= 69069;
    t *= 69069;
    return .5 * (double)t / (double)(unsigned int)0x80000000;
}

static double
trand5(i)
    int i;
{
    unsigned int t = 2 * i + 1;
    t *= 1664525;
    t *= 1664525;
    t *= 1664525;
    return .5 * (double)t / (double)(unsigned int)0x80000000;
}


/***************  sample  ***************/

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

	/**  right-hand side  **/

	*b = trand3(i);

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


/***************  honsen1  ***************/

/*
 *  2-dimensional Poisson Equation with Dirichlet Condition
 *
 *  Spiral Ordering
 *
 *  Assume SIZE_X == SIZE_Y, and SIZE_X is odd.
 */

#define SIZE_X		1413
#define SIZE_Y		(SIZE_X)
#define MAX_NZERO	5
#define RESMAX		1e-5

static void
r_to_ij(n, r, i, j)
    int n, r;
    int *i, *j;
{
    int c, d, e, f;
    double d_n = (double)n;

    c = (int)floor(.5 * (d_n - sqrt(d_n * d_n - r)));
    d = r - 4 * c * (n - c);
    e = n - 2 * c - 1;
    if (e > 0)
	f = d / e;
    else
	f = 0;

    switch (f) {
	char errmsg[255];
    case 0:  /**  upper edge  **/
	*i = c;
	*j = c + d;
	break;
    case 1:  /**  right edge  **/
	*i = c + d - e;
	*j = n - 1 - c;
	break;
    case 2:  /**  lower edge  **/
	*i = n - 1 - c;
	*j = n - 1 - c - (d - 2 * e);
	break;
    case 3:  /**  left edge  **/
	*i = n - 1 - c - (d - 3 * e);
	*j = c;
	break;
    default:
	sprintf(errmsg, "r_to_ij: r = %d, f = %d", r, f);
	perror(errmsg);
	exit(1);
    }
}

#define MAX(i, j)  ((i) > (j) ? (i) : (j))

static int
ij_to_r(n, i, j)
    int n, i, j;
{
    int i1, j1, c;
    int n2 = n / 2;

    i1 = i - n2;
    j1 = j - n2;
    c = n2 - MAX(abs(i1), abs(j1));

    if (i - j <= 0) {
	return 4 * c * (n - c) + i + j - 2 * c;
    }
    else {
	return 4 * (c + 1)* (n - c - 1) - (i + j - 2 * c);
    }
}

#undef MAX

static int
honsen1(i, ja, a, rhs)
    int i, *ja;
    double *a, *rhs;
{
    double inv_h = (SIZE_X + 1) * (SIZE_X + 1);
    double d0 = 4.0 * inv_h;
    double d1 = -1.0 * inv_h;
    int N = SIZE_X * SIZE_Y;

    if (i == -1) {
	ja[0] = N;
	ja[1] = N + 2 * (N - SIZE_Y)
	    + 2 * (N - SIZE_X);
	ja[2] = MAX_NZERO;
	*rhs  = RESMAX;

	assert(N <= MATMAX);

	return 0;
    }
    else if (i >= 0 && i < N) {
	double b[MAX_NZERO];
	int jb[MAX_NZERO];
	int perm[MAX_NZERO];

	int k = 0;
	int s, t;
	int j;

	r_to_ij(SIZE_X, i, &s, &t);

	/**  Diagonal elements  **/

	b[k] = d0;
#if IDX_BEGINS_WITH_ZERO
	jb[k] = i;
#else
	jb[k] = i + 1;
#endif
	++k;

	/**  Subdiagonal elements  **/

	if (s - 1 >= 0) {
	    b[k] = d1;
#if IDX_BEGINS_WITH_ZERO
	    jb[k] = ij_to_r(SIZE_X, s - 1, t);
#else
	    jb[k] = ij_to_r(SIZE_X, s - 1, t) + 1;
#endif
	    ++k;
	}
	if (t - 1 >= 0) {
	    b[k] = d1;
#if IDX_BEGINS_WITH_ZERO
	    jb[k] = ij_to_r(SIZE_X, s, t - 1);
#else
	    jb[k] = ij_to_r(SIZE_X, s, t - 1) + 1;
#endif
	    ++k;
	}
	if (s + 1 < SIZE_X) {
	    b[k] = d1;
#if IDX_BEGINS_WITH_ZERO
	    jb[k] = ij_to_r(SIZE_X, s + 1, t);
#else
	    jb[k] = ij_to_r(SIZE_X, s + 1, t) + 1;
#endif
	    ++k;
	}
	if (t + 1 < SIZE_Y) {
	    b[k] = d1;
#if IDX_BEGINS_WITH_ZERO
	    jb[k] = ij_to_r(SIZE_X, s, t + 1);
#else
	    jb[k] = ij_to_r(SIZE_X, s, t + 1) + 1;
#endif
	    ++k;
	}

	sort_PSC98(k, jb, perm);

	for (j = 0; j < k; ++j) {
	    ja[j] = jb[perm[j]];
	    a[j] = b[perm[j]];
	}

	for (; k < NZMAX; ++k) {
	    ja[k] = -1;
	    a[k] = 0.0;
	}

	/**  right-hand side  **/

	*rhs = trand1(i);

	return 0;
    }
    else {
	perror("genmat");

	return 1;
    }
}

#undef SIZE_X
#undef SIZE_Y
#undef MAX_NZERO
#undef RESMAX


/***************  honsen2  ***************/

/*
 *  2-dimensional Poisson Equation with Dirichlet Condition
 *
 *  Spiral Ordering
 *
 *  Assume SIZE_X == SIZE_Y, and SIZE_X is odd.
 */

#define SIZE_X		511
#define SIZE_Y		(SIZE_X)
#define MAX_NZERO	5
#define RESMAX		1e-5

/*
 *  +----------+
 *  | k1       |
 *  |    +--+  |
 *  |    |k2|  |
 *  |    +--+  |
 *  |          |
 *  +----------+
 */

static void
get_coeff2(i, c)
    int i;
    double *c;
{
    double k1 = 1;
    double k2 = 100;
    int ix, iy;
    int cx1, cy1, cx2, cy2;

    iy = i / SIZE_X;
    ix = i % SIZE_X;

    cx1 = SIZE_X / 2;
    cy1 = SIZE_Y / 2;
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

/*
 *  Spiral problem
 */

static int
honsen2(i, ja, a, rhs)
    int i, *ja;
    double *a, *rhs;
{
    double inv_h = (SIZE_X + 1) * (SIZE_X + 1);
    int N = SIZE_X * SIZE_Y;

    if (i == -1) {
	ja[0] = N;
	ja[1] = N + 2 * (N - SIZE_Y)
	    + 2 * (N - SIZE_X);
	ja[2] = MAX_NZERO;
	*rhs  = RESMAX;

	assert(N <= MATMAX);

	return 0;
    }
    else if (i >= 0 && i < N) {
	double b[MAX_NZERO];
	int jb[MAX_NZERO];
	int perm[MAX_NZERO];

	int k = 0;
	int s, t;
	int j;
	double c[4];
	double c0, c1, c2, c3;

	r_to_ij(SIZE_X, i, &s, &t);

	/*
	 *  c[2] | c[3]
	 *  -----+-----
	 *  c[0] | c[1]
	 */

	get_coeff2(t * SIZE_X + s, c);

	c0 = c[0] * inv_h;
	c1 = c[1] * inv_h;
	c2 = c[2] * inv_h;
	c3 = c[3] * inv_h;

	/**  Diagonal elements  **/

	b[k] = c0 + c1 + c2 + c3;
#if IDX_BEGINS_WITH_ZERO
	jb[k] = i;
#else
	jb[k] = i + 1;
#endif
	++k;

	/**  Subdiagonal elements  **/

	if (s - 1 >= 0) {
	    b[k] = -.5 * (c0 + c2);
#if IDX_BEGINS_WITH_ZERO
	    jb[k] = ij_to_r(SIZE_X, s - 1, t);
#else
	    jb[k] = ij_to_r(SIZE_X, s - 1, t) + 1;
#endif
	    ++k;
	}
	if (t - 1 >= 0) {
	    b[k] = -.5 * (c0 + c1);
#if IDX_BEGINS_WITH_ZERO
	    jb[k] = ij_to_r(SIZE_X, s, t - 1);
#else
	    jb[k] = ij_to_r(SIZE_X, s, t - 1) + 1;
#endif
	    ++k;
	}
	if (s + 1 < SIZE_X) {
	    b[k] = -.5 * (c1 + c3);
#if IDX_BEGINS_WITH_ZERO
	    jb[k] = ij_to_r(SIZE_X, s + 1, t);
#else
	    jb[k] = ij_to_r(SIZE_X, s + 1, t) + 1;
#endif
	    ++k;
	}
	if (t + 1 < SIZE_Y) {
	    b[k] = -.5 * (c2 + c3);
#if IDX_BEGINS_WITH_ZERO
	    jb[k] = ij_to_r(SIZE_X, s, t + 1);
#else
	    jb[k] = ij_to_r(SIZE_X, s, t + 1) + 1;
#endif
	    ++k;
	}

	sort_PSC98(k, jb, perm);

	for (j = 0; j < k; ++j) {
	    ja[j] = jb[perm[j]];
	    a[j] = b[perm[j]];
	}

	for (; k < NZMAX; ++k) {
	    ja[k] = -1;
	    a[k] = 0.0;
	}

	/**  right-hand side  **/

	*rhs = trand2(i);

	return 0;
    }
    else {
	perror("genmat");

	return 1;
    }
}

#undef SIZE_X
#undef SIZE_Y
#undef MAX_NZERO
#undef RESMAX


/***************  honsen3  ***************/

/*
 *  2-dimensional Poisson Equation with Dirichlet Condition
 *
 *  9-point stencil
 */

#define SIZE_X		1413
#define SIZE_Y		1413
#define PROB_SIZE	(SIZE_X * SIZE_Y)
#define MAX_NZERO	9
#define RESMAX		1e-5

static double
rhs3(i)
    int i;
{
    double q1, q2, q3, q4, q5;
    if (i - SIZE_X >= 0)
	q1 = trand3(i - SIZE_X);
    else
	q1 = 0.0;

    if (i % SIZE_X > 0)
	q2 = trand3(i - 1);
    else
	q2 = 0.0;

    q3 = trand3(i);

    if ((i + 1) % SIZE_X > 0)
	q4 = trand3(i + 1);
    else
	q4 = 0.0;

    if (i + SIZE_X < PROB_SIZE)
	q5 = trand3(i + SIZE_X);
    else
	q5 = 0.0;

    return q1 + q2 + 8.0 * q3 + q4 + q5;
}

static int
honsen3(i, ja, a, b)
    int i, *ja;
    double *a, *b;
{
    double inv_h = (SIZE_X + 1) * (SIZE_X + 1);
    double d0 = 40.0 * inv_h;
    double d1 = -8.0 * inv_h;
    double d2 = -2.0 * inv_h;

    if (i == -1) {
	ja[0] = PROB_SIZE;
	ja[1] = PROB_SIZE
	    + 2 * (PROB_SIZE - SIZE_Y)
	    + 2 * (PROB_SIZE - SIZE_X)
	    + 4 * (PROB_SIZE - SIZE_X - SIZE_Y + 1);
	ja[2] = MAX_NZERO;
	*b    = RESMAX;

	assert(PROB_SIZE <= MATMAX);

	return 0;
    }
    else if (i >= 0 && i < PROB_SIZE) {
	int k = 0;

	if (i - SIZE_X >= 0 && i % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - SIZE_X - 1;
#else
	    ja[k] = i - SIZE_X - 1 + 1;
#endif
	    a[k] = d2;
	    ++k;
	}
	if (i - SIZE_X >= 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - SIZE_X;
#else
	    ja[k] = i - SIZE_X + 1;
#endif
	    a[k] = d1;
	    ++k;
	}
	if (i - SIZE_X >= 0 && (i + 1) % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - SIZE_X + 1;
#else
	    ja[k] = i - SIZE_X + 1 + 1;
#endif
	    a[k] = d2;
	    ++k;
	}
	if (i % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - 1;
#else
	    ja[k] = i - 1 + 1;
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
	    ja[k] = i + 1 + 1;
#endif
	    a[k] = d1;
	    ++k;
	}
	if (i + SIZE_X < PROB_SIZE && i % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + SIZE_X - 1;
#else
	    ja[k] = i + SIZE_X - 1 + 1;
#endif
	    a[k] = d2;
	    ++k;
	}
	if (i + SIZE_X < PROB_SIZE) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + SIZE_X;
#else
	    ja[k] = i + SIZE_X + 1;
#endif
	    a[k] = d1;
	    ++k;
	}
	if (i + SIZE_X < PROB_SIZE && (i + 1) % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + SIZE_X + 1;
#else
	    ja[k] = i + SIZE_X + 1 + 1;
#endif
	    a[k] = d2;
	    ++k;
	}
	for (; k < NZMAX; ++k) {
	    ja[k] = -1;
	    a[k] = 0.0;
	}

	/**  right-hand side  **/

	*b = rhs3(i);

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


/***************  honsen4  ***************/

/*
 *  3-dimensional Poisson Equation with Dirichlet Condition
 *
 *  Non-uniform diffusion constant
 */

#define SIZE_X		101
#define SIZE_Y		101
#define SIZE_Z		101
#define PROB_SIZE	(SIZE_X * SIZE_Y * SIZE_Z)
#define SIZE_XY		(SIZE_X * SIZE_Y)
#define MAX_NZERO	7
#define	RESMAX		1e-5

/*
 *  +----------+
 *  | k1       |
 *  |    +--+  |
 *  |    |k2|  |
 *  |    +--+  |
 *  |          |
 *  +----------+
 */

static void
get_coeff4(i, c)
    int i;
    double *c;
{
    double k1 = 1;
    double k2 = 100;
    int ix, iy, iz;
    int cx1, cy1, cz1, cx2, cy2, cz2;

    iz = i / SIZE_XY;
    iy = (i % SIZE_XY) / SIZE_X;
    ix = i % SIZE_XY % SIZE_X;

    cx1 = SIZE_X / 2;
    cy1 = SIZE_Y / 2;
    cz1 = SIZE_Z / 2;
    cx2 = 3 * SIZE_X / 4;
    cy2 = 3 * SIZE_Y / 4;
    cz2 = 3 * SIZE_Z / 4;

    /*
     *  c[2] | c[3]
     *  -----+-----
     *  c[0] | c[1]
     *
     *     c[6] | c[7]
     *     -----+-----
     *     c[4] | c[5]
     */

    if ((ix > cx1 && iy > cy1 && iz > cz1)
	&& (ix <= cx2 && iy <= cy2 && iz <= cz2))
	c[0] = k2;
    else
	c[0] = k1;
	
    if ((ix >= cx1 && iy > cy1 && iz > cz1)
	&& (ix < cx2 && iy <= cy2 && iz <= cz2))
	c[1] = k2;
    else
	c[1] = k1;

    if ((ix > cx1 && iy >= cy1 && iz > cz1)
	&& (ix <= cx2 && iy < cy2 && iz <= cz2))
	c[2] = k2;
    else
	c[2] = k1;

    if ((ix >= cx1 && iy >= cy1 && iz > cz1)
	&& (ix < cx2 && iy < cy2 && iz <= cz2))
	c[3] = k2;
    else
	c[3] = k1;

    if ((ix > cx1 && iy > cy1 && iz >= cz1)
	&& (ix <= cx2 && iy <= cy2 && iz < cz2))
	c[4] = k2;
    else
	c[4] = k1;
	
    if ((ix >= cx1 && iy > cy1 && iz >= cz1)
	&& (ix < cx2 && iy <= cy2 && iz < cz2))
	c[5] = k2;
    else
	c[5] = k1;

    if ((ix > cx1 && iy >= cy1 && iz >= cz1)
	&& (ix <= cx2 && iy < cy2 && iz < cz2))
	c[6] = k2;
    else
	c[6] = k1;

    if ((ix >= cx1 && iy >= cy1 && iz >= cz1)
	&& (ix < cx2 && iy < cy2 && iz < cz2))
	c[7] = k2;
    else
	c[7] = k1;
}

static int
honsen4(i, ja, a, b)
    int i, *ja;
    double *a, *b;
{
    double inv_h = (SIZE_X + 1) * (SIZE_X + 1);

    if (i == -1) {
	ja[0] = PROB_SIZE;
	ja[1] = PROB_SIZE
	    + 2 * (PROB_SIZE - SIZE_XY)
	    + 2 * (PROB_SIZE - SIZE_Y * SIZE_Z)
	    + 2 * (PROB_SIZE - SIZE_Z * SIZE_X);
	ja[2] = MAX_NZERO;
	*b    = RESMAX;

	assert(PROB_SIZE <= MATMAX);

	return 0;
    }
    else if (i >= 0 && i < PROB_SIZE) {
	int k = 0;
	double c[8];
	double c0, c1, c2, c3;
	double c4, c5, c6, c7;

	/*
	 *  c[2] | c[3]
	 *  -----+-----
	 *  c[0] | c[1]
	 *
	 *     c[6] | c[7]
	 *     -----+-----
	 *     c[4] | c[5]
	 */

	get_coeff4(i, c);

	c0 = c[0] * inv_h;
	c1 = c[1] * inv_h;
	c2 = c[2] * inv_h;
	c3 = c[3] * inv_h;
	c4 = c[4] * inv_h;
	c5 = c[5] * inv_h;
	c6 = c[6] * inv_h;
	c7 = c[7] * inv_h;

	if (i >= SIZE_XY) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - SIZE_XY;
#else
	    ja[k] = i - SIZE_XY + 1;
#endif
	    a[k] = -.25 * (c0 + c1 + c2 + c3);
	    ++k;
	}
	if (i % SIZE_XY >= SIZE_X) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - SIZE_X;
#else
	    ja[k] = i - SIZE_X + 1;
#endif
	    a[k] = -.25 * (c0 + c1 + c4 + c5);
	    ++k;
	}
	if (i % SIZE_XY % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - 1;
#else
	    ja[k] = i - 1 + 1;
#endif
	    a[k] = -.25 * (c0 + c2 + c4 + c6);
	    ++k;
	}
#if IDX_BEGINS_WITH_ZERO
	ja[k] = i;
#else
	ja[k] = i + 1;
#endif
	a[k] = c0 + c1 + c2 + c3 + c4 + c5 + c6 + c7;
	++k;
	if ((i + 1) % SIZE_XY % SIZE_X > 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + 1;
#else
	    ja[k] = i + 1 + 1;
#endif
	    a[k] = -.25 * (c1 + c3 + c5 + c7);
	    ++k;
	}
	if (i % SIZE_XY < SIZE_XY - SIZE_X) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + SIZE_X;
#else
	    ja[k] = i + SIZE_X + 1;
#endif
	    a[k] = -.25 * (c2 + c3 + c6 + c7);
	    ++k;
	}
	if (i < PROB_SIZE - SIZE_XY) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + SIZE_XY;
#else
	    ja[k] = i + SIZE_XY + 1;
#endif
	    a[k] = -.25 * (c4 + c5 + c6 + c7);
	    ++k;
	}
	for (; k < NZMAX; ++k) {
	    ja[k] = -1;
	    a[k] = 0.0;
	}

	/**  right-hand side  **/

	*b = trand4(i);

	return 0;
    }
    else {
	perror("genmat");

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


/***************  honsen5  ***************/

/*
 *  1-dimensional Poisson Equation with Dirichlet Condition
 *
 *  Non-uniform diffusion constant
 */

#define SIZE_X	        24575
#define PROB_SIZE	(SIZE_X)
#define MAX_NZERO	3
#define	RESMAX		1e-3

/*
 *     k1  k2 k1
 *  +-----+--+--+
 */

static void
get_coeff5(i, c)
    int i;
    double *c;
{
    double k1 = 1;
    double k2 = 10;
    int cx1, cx2;

    cx1 = SIZE_X / 2;
    cx2 = 3 * SIZE_X / 4;

    if ((i > cx1) && (i <= cx2))
	c[0] = k2;
    else
	c[0] = k1;
	
    if ((i >= cx1) && (i < cx2))
	c[1] = k2;
    else
	c[1] = k1;
}

int
honsen5(i, ja, a, b)
    int i, *ja;
    double *a, *b;
{
    double inv_h = (SIZE_X + 1) * (SIZE_X + 1);

    if (i == -1) {
	ja[0] = PROB_SIZE;
	ja[1] = PROB_SIZE + 2 * (PROB_SIZE - 1);
	ja[2] = MAX_NZERO;
	*b    = RESMAX;

	assert(PROB_SIZE <= MATMAX);

	return 0;
    }
    else if (i >= 0 && i < PROB_SIZE) {
	int k = 0;
	double c[2];
	double c0, c1;

	/*
	 *  c[0] | c[1]
	 */

	get_coeff5(i, c);

	c0 = c[0] * inv_h;
	c1 = c[1] * inv_h;

	if (i - 1 >= 0) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i - 1;
#else
	    ja[k] = i - 1 + 1;
#endif
	    a[k] = -c0;
	    ++k;
	}
#if IDX_BEGINS_WITH_ZERO
	ja[k] = i;
#else
	ja[k] = i + 1;
#endif
	a[k] = c0 + c1;
	++k;
	if (i + 1 < PROB_SIZE) {
#if IDX_BEGINS_WITH_ZERO
	    ja[k] = i + 1;
#else
	    ja[k] = i + 1 + 1;
#endif
	    a[k] = -c1;
	    ++k;
	}
	for (; k < NZMAX; ++k) {
	    ja[k] = -1;
	    a[k] = 0.0;
	}

	/**  right-hand side  **/

	*b = trand5(i);

	return 0;
    }
    else {
	perror("genmat");

	return 1;
    }
}

#undef SIZE_X
#undef PROB_SIZE
#undef MAX_NZERO
#undef RESMAX

