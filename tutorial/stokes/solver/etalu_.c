/* etalu_.f -- translated by f2c (version 19950920).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;


/* Subroutine */ int etalu_(a, ia, l, l2, b, nwk, m, eps, num, x, e, ie, we, 
	ne, ip, is, ier)
doublereal *a;
integer *ia, *l, *l2;
doublereal *b;
integer *nwk, *m;
doublereal *eps;
integer *num;
doublereal *x, *e;
integer *ie;
doublereal *we;
integer *ne, *ip, *is, *ier;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static integer irow;
    static doublereal work;
    static integer i__, k;
    extern /* Subroutine */ int updtl_(), updtu_();
    static integer jk;
    static doublereal val, piv, eps1;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___9 = { 0, 6, 0, 0, 0 };


/* ********************************************************************** 
*/
/*  GAUSS METHOD FOR A NON-SYMMETRIC SPARSE MATRIX.                    * 
*/
/*                                                                     * 
*/
/*  PARAMETERS:                                                        * 
*/
/*   ON ENTRY:                                                         * 
*/
/*     A      THE ARRAY WHICH CONTAINS NON-ZERO ELEMENTS IN COLUMN-WISE* 
*/
/*     IA     THE ARRAY WHICH HAS CORRESPONDING ROW INDEX WITH ARRAY A.* 
*/
/*     L      THE LEADING DIMENSION OF THE ARRAY A.                    * 
*/
/*     L2     THE LEADING DIMENSION OF THE ARRAY E.                    * 
*/
/*     B      THE RIGHT HAND SIDE VECTOR.                              * 
*/
/*     NWK    ACCUMULATED SUM OF NON-ZERO ELEMENTS IN EACH COLUMN OF A.* 
*/
/*     M      THE ORDER OF THE MATRIX A.                               * 
*/
/*     EPS    THE TOLERANCE FOR PIVOTAL ELEMENTS.                      * 
*/
/*   ON RETURN:                                                        * 
*/
/*     X      THE SOLUTION VECTOR.                                     * 
*/
/*     E      THE ARRAY WHICH CONTAINS NON ZERO ELEMENTS IN ETA-VECTORS* 
*/
/*            IN COLUMN-WISE.                                          * 
*/
/*     IE     THE ARRAY WHICH HAS CORRESPONDING ROW INDEX WITH ARRAY E.* 
*/
/*     NETA   ACCUMULATED SUM OF NON-ZERO ELEMENTS IN EACH COLUMN OF E.* 
*/
/*            EACH COLUMN HAS TWO ENTRIES FOR U-ETA AND L-ETA.         * 
*/
/*     IP     THE ARRAY WHICH HAS PIVOTAL ROWS.                        * 
*/
/*     IS     THE NUMBER OF DEGENERATED COLUMNS.                       * 
*/
/*     IER    THE ERROR CODE. IF IER=0, NORMAL RETURN.                 * 
*/
/*   OTHERS:  WORKING PARAMETERS.                                      * 
*/
/*                                                                     * 
*/
/*  COPYRIFHT:    TSUTOMU OGUNI      SEP. 1 1992      VER. 2           * 
*/
/* ********************************************************************** 
*/

    /* Parameter adjustments */
    --ia;
    --a;
    --ie;
    --e;
    --ip;
    ++ne;
    --we;
    --x;
    --b;

    /* Function Body */
    if (*l < *m) {
	s_wsle(&io___1);
	do_lio(&c__9, &c__1, "(SUBR. ETALU) INVALID ARGUMENT. ", 32L);
	do_lio(&c__3, &c__1, (char *)&(*l), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&(*m), (ftnlen)sizeof(integer));
	e_wsle();
	*ier = 2;
	return 0;
    }
    eps1 = *eps * .01;
    *ier = 0;
    *is = 0;
    ne[-1] = 0;
    ne[0] = 0;
    *num = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	ip[i__] = 0;
    }

    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L20: */
	    we[i__] = 0.;
	}
	i__2 = nwk[k];
	for (i__ = nwk[k - 1] + 1; i__ <= i__2; ++i__) {
/* L30: */
	    we[ia[i__]] = a[i__];
	}
	if (k != 1) {
	    jk = k - 1;
	    updtl_(&e[1], &ie[1], &ne[-1], &we[1], &ip[1], l2, &jk, m);
/*  COMPUTATION OF ETA */
	}
/*  SELECTION OF PIVOT */
	irow = 0;
	val = 0.;
	i__2 = *m;
	for (i__ = k; i__ <= i__2; ++i__) {
	    work = (d__1 = we[i__], abs(d__1));
	    if (work > val) {
		val = work;
		irow = i__;
	    }
/* L60: */
	}
	if (val <= *eps) {
	    s_wsle(&io___9);
	    do_lio(&c__9, &c__1, "(SUBR. ETALU) STOP WITH ZERO PIVOT", 34L);
	    do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsle();
	    *ier = 1;
	    return 0;
	}
	ip[k] = irow;
	if (irow != k) {
	    work = we[k];
	    we[k] = we[irow];
	    we[irow] = work;
	}
/*      WRITE(*,*) (WE(I),I=1,M) */
/*  GENERATION OF ETA */
	piv = -1. / we[k];
	we[k] = -1.;
	i__2 = *m;
	for (i__ = k; i__ <= i__2; ++i__) {
	    we[i__] *= piv;
	    if ((d__1 = we[i__], abs(d__1)) >= eps1) {
		++(*num);
		e[*num] = we[i__];
		ie[*num] = i__;
	    }
/* L70: */
	}
	ne[(k << 1) - 1] = *num;
	if (k == 1) {
	    ne[2] = ne[1];
	} else {
	    i__2 = k - 1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if ((d__1 = we[i__], abs(d__1)) >= eps1) {
		    ++(*num);
		    e[*num] = -we[i__];
		    ie[*num] = i__;
		}
/* L80: */
	    }
	    ne[k * 2] = *num;
	}

/* L90: */
    }
/*  COMPUTATION OF X */
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
/* L100: */
	x[k] = b[k];
    }
    updtl_(&e[1], &ie[1], &ne[-1], &x[1], &ip[1], l2, m, m);
    updtu_(&e[1], &ie[1], &ne[-1], &x[1], l2, m);

    return 0;
} /* etalu_ */


/* Subroutine */ int updtl_(e, ie, ne, we, ip, l2, k, m)
doublereal *e;
integer *ie, *ne;
doublereal *we;
integer *ip, *l2, *k, *m;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer irow;
    static doublereal work;
    static integer i__, j;


    /* Parameter adjustments */
    --ie;
    --e;
    --ip;
    --we;
    ++ne;

    /* Function Body */
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	irow = ip[j];
	if (irow != j) {
	    work = we[j];
	    we[j] = we[irow];
	    we[irow] = work;
	}
	work = we[j];
	we[j] = 0.;
	i__2 = ne[(j << 1) - 1];
	for (i__ = ne[(j << 1) - 2] + 1; i__ <= i__2; ++i__) {
/* L10: */
	    we[ie[i__]] += e[i__] * work;
	}
/* L20: */
    }
    return 0;
} /* updtl_ */


/* Subroutine */ int updtu_(e, ie, ne, we, l2, m)
doublereal *e;
integer *ie, *ne;
doublereal *we;
integer *l2, *m;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal work;
    static integer i__, j;


    /* Parameter adjustments */
    --ie;
    --e;
    --we;
    ++ne;

    /* Function Body */
    for (j = *m; j >= 2; --j) {
	if (ne[(j << 1) - 1] != ne[j * 2]) {
	    work = we[j];
	    i__1 = ne[j * 2];
	    for (i__ = ne[(j << 1) - 1] + 1; i__ <= i__1; ++i__) {
/* L50: */
		we[ie[i__]] += e[i__] * work;
	    }
	}
/* L40: */
    }
    return 0;
/*  END OF ETALU */
} /* updtu_ */

