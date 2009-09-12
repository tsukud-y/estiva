/* pcg_.f -- translated by f2c (version 19950920).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;


/* Subroutine */ int pcg_(d__, a, ia, n, n1, nl, b, eps, itr, s, x, dd, p, q, 
	r__, m, ier)
doublereal *d__, *a;
integer *ia, *n, *n1, *nl;
doublereal *b, *eps;
integer *itr;
doublereal *s, *x, *dd, *p, *q, *r__;
integer *m, *ier;
{
    /* System generated locals */
    integer a_dim1, a_offset, ia_dim1, ia_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();
    double sqrt();

    /* Local variables */
    static doublereal beta;
    static integer i__, j, k, l;
    static doublereal alpha, w, y, c1, c2, c3, x1, x2, th;
    static integer nn;
    static doublereal ss, res;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };


/* ********************************************************************* 
*/
/*  ICCG AND MICCG METHOD FOR FINITE ELEMENT METHOD.                  * */
/*                                                                    * */
/*   PARAMETERS:                                                      * */
/*    ON ENTRY:                                                       * */
/*      D      1-DIM. ARRAY CONTAINING DIAGONAL ELEMENTS OF MATRIX A. * */
/*      A      2-DIM. ARRAY CONTAINING NON-AERO ELEMENTS OF LOWER     * */
/*             PART OF THE MATRIX A EXCEPT DIAGONALS.                 * */
/*      IA     2-DIM. ARRAY CONTAINING THE COLUMN INDEX OF ELEMENTS   * */
/*             IN TH ARRAY A.                                         * */
/*      N      THE ORDER OF THE MATRIX A.                             * */
/*      N1     THE LEADING DIMENSION OF THE ARRAY A.                  * */
/*      NL     MORE THAN MAXIMUM NUMBER OF NON-ZERO ELEMENTS IN EACH  * */
/*             ROW OF THE ARRAY A.                                    * */
/*      B      1-DIM. ARRAY CONTAINING RIGHT HAND SIDE VECTOR.        * */
/*      EPS    THE TOLERANCE FOR CONVERGENCE.                         * */
/*      ITR    MAXIMUM NUMBER OF ITERATIONS.                          * */
/*      S      THE SIGMA VALUE WHICH SPECIFIES A METHOD TO BE USED.   * */
/*   ON RETURN:                                                       * */
/*      X      1-DIM. ARRAY CONTAINING THE SOLUTION VECTOR.           * */
/*      EPS    RELATIVE ERROR AT RETURN.                              * */
/*      ITR    NUMBER OF ITERATIONS AT RETURN.                        * */
/*      IER    THE ERROR CODE. IF IER=0, NORMAL RETURN.               * */
/*   OTHERS:   WORKING PARAMETERS:                                    * */
/*                                                                    * */
/*  COPYRIGHT:   TSUTOMU OGUNI    FEB. 1 1993   VER. 2                * */
/* ********************************************************************* 
*/

    /* Parameter adjustments */
    --m;
    --r__;
    --b;
    --d__;
    ia_dim1 = *n1;
    ia_offset = ia_dim1 + 1;
    ia -= ia_offset;
    a_dim1 = *n1;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    *ier = 0;
    if (*n1 < *n || *s < (float)0.) {
	s_wsle(&io___1);
	do_lio(&c__9, &c__1, "(SUBR. PCG) INVALID ARGUMENT. ", 30L);
	do_lio(&c__3, &c__1, (char *)&(*n1), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&(*nl), (ftnlen)sizeof(integer));
	do_lio(&c__5, &c__1, (char *)&(*s), (ftnlen)sizeof(doublereal));
	e_wsle();
	*ier = 2;
	return 0;
    }
/*  INITIALIZATION */
    th = 1.;
    if (*s > (float)0. && *s < (float)1.) {
	th = *s;
	*s = 1.;
    }
    i__1 = *n << 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L5: */
	m[i__] = 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dd[i__] = 0.;
	i__2 = *nl;
	for (j = 1; j <= i__2; ++j) {
	    k = ia[i__ + j * ia_dim1];
	    if (k != 0) {
		++m[*n + k];
		ia[k + (m[*n + k] + *nl) * ia_dim1] = i__;
		a[k + (m[*n + k] + *nl) * a_dim1] = a[i__ + j * a_dim1];
		++m[i__];
	    }
/* L10: */
	}
    }
    x[0] = 0.;
    dd[0] = 0.;
    p[0] = 0.;
    q[0] = 0.;
    if (*s != (float)0.) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    w = d__[i__] * *s;
	    i__1 = m[i__];
	    for (k = 1; k <= i__1; ++k) {
		nn = ia[i__ + k * ia_dim1];
		ss = a[i__ + k * a_dim1];
		if (nn != 0) {
		    i__3 = *nl + m[nn + *n];
		    for (j = *nl + 1; j <= i__3; ++j) {
			if (ia[nn + j * ia_dim1] != i__) {
			    ss += a[nn + j * a_dim1] * th;
			}
/* L16: */
		    }
		}
/* L14: */
		w -= a[i__ + k * a_dim1] * ss * dd[nn];
	    }
/* L12: */
	    dd[i__] = 1. / w;
	}
    } else {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ss = d__[i__];
	    i__1 = m[i__];
	    for (k = 1; k <= i__1; ++k) {
/* L20: */
/* Computing 2nd power */
		d__1 = a[i__ + k * a_dim1];
		ss -= d__1 * d__1 * dd[ia[i__ + k * ia_dim1]];
	    }
/* L30: */
	    dd[i__] = 1. / ss;
	}
    }
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__] = d__[i__] * x[i__];
	i__1 = m[i__];
	for (j = 1; j <= i__1; ++j) {
/* L40: */
	    q[i__] += a[i__ + j * a_dim1] * x[ia[i__ + j * ia_dim1]];
	}
	i__1 = *nl + m[i__ + *n];
	for (j = *nl + 1; j <= i__1; ++j) {
/* L45: */
	    q[i__] += a[i__ + j * a_dim1] * x[ia[i__ + j * ia_dim1]];
	}
/* L50: */
    }
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L60: */
	r__[i__] = b[i__] - r__[i__];
    }
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	p[i__] = r__[i__];
	i__1 = m[i__];
	for (k = 1; k <= i__1; ++k) {
/* L70: */
	    p[i__] -= a[i__ + k * a_dim1] * p[ia[i__ + k * ia_dim1]];
	}
/* L80: */
	p[i__] = dd[i__] * p[i__];
    }
    for (i__ = *n; i__ >= 1; --i__) {
	ss = 0.;
	i__2 = *nl + m[i__ + *n];
	for (k = *nl + 1; k <= i__2; ++k) {
/* L90: */
	    ss += a[i__ + k * a_dim1] * p[ia[i__ + k * ia_dim1]];
	}
/* L100: */
	p[i__] -= dd[i__] * ss;
    }
    c1 = 0.;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L110: */
	c1 += r__[i__] * p[i__];
    }
/*  ITERATION PHASE */
    i__2 = *itr;
    for (l = 1; l <= i__2; ++l) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    q[i__] = d__[i__] * p[i__];
	    i__3 = m[i__];
	    for (j = 1; j <= i__3; ++j) {
/* L130: */
		q[i__] += a[i__ + j * a_dim1] * p[ia[i__ + j * ia_dim1]];
	    }
	    i__3 = *nl + m[i__ + *n];
	    for (j = *nl + 1; j <= i__3; ++j) {
/* L135: */
		q[i__] += a[i__ + j * a_dim1] * p[ia[i__ + j * ia_dim1]];
	    }
/* L140: */
	}
	c2 = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L150: */
	    c2 += p[i__] * q[i__];
	}
	if (c2 == (float)0.) {
	    *ier = 3;
	    *itr = l;
	    goto L300;
	}
	alpha = c1 / c2;
	x1 = 0.;
	x2 = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    y = x[i__];
	    x[i__] += alpha * p[i__];
	    r__[i__] -= alpha * q[i__];
	    x1 += y * y;
/* L160: */
/* Computing 2nd power */
	    d__1 = x[i__] - y;
	    x2 += d__1 * d__1;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    q[i__] = r__[i__];
	    i__3 = m[i__];
	    for (k = 1; k <= i__3; ++k) {
/* L180: */
		q[i__] -= a[i__ + k * a_dim1] * q[ia[i__ + k * ia_dim1]];
	    }
/* L190: */
	    q[i__] = dd[i__] * q[i__];
	}
	for (i__ = *n; i__ >= 1; --i__) {
	    i__1 = *nl + m[i__ + *n];
	    for (k = *nl + 1; k <= i__1; ++k) {
/* L220: */
		q[i__] -= dd[i__] * a[i__ + k * a_dim1] * q[ia[i__ + k * 
			ia_dim1]];
	    }
	}
	if (x1 != (float)0.) {
	    res = sqrt(x2 / x1);
	    if (res < *eps) {
		*itr = l;
		*eps = res;
		*ier = 0;
		if (th != (float)1.) {
		    *s = th;
		}
		return 0;
	    }
	}
	c3 = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L230: */
	    c3 += r__[i__] * q[i__];
	}
	if (c1 == (float)0.) {
	    *itr = l;
	    *ier = 4;
	    goto L300;
	}
	beta = c3 / c1;
	c1 = c3;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    p[i__] = q[i__] + beta * p[i__];
/* L240: */
	}

/* L250: */
    }
    *ier = 1;
    s_wsle(&io___19);
    do_lio(&c__9, &c__1, "(SUBR. PCG) NO CONVERGENCE.", 27L);
    e_wsle();
L300:
    *eps = res;
    if (th != (float)1.) {
	*s = th;
    }
    return 0;
/*  END OF PCG */
} /* pcg_ */

