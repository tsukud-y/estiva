/* STWART.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <stdio.h>
#include "f2c.h"

/* Table of constant values */

static long c__9 = 9;
static long c__1 = 1;
static long c__3 = 3;


/* Subroutine */ int stwart_(ia, l, ma, m, r__, c__, ir, ic, jrow, jcol, ip, 
	jp, kerns, mend, iw, lg, ier)
long *ia, *l, *ma, *m, *r__, *c__, *ir, *ic, *jrow, *jcol, *ip, *jp, *
	kerns, *mend, *iw, *lg, *ier;
{
    /* System generated locals */
    long i__1, i__2, i__3;

    /* Builtin functions */
    long s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static long kern, irow, i__, j, k, n, mm, is;
    extern /* Subroutine */ int single_(), forful_();
    static long min__;

    /* Fortran I/O blocks */
    static cilist io___8 = { 0, 6, 0, 0, 0 };


/* ********************************************************************** */
/*  STWART METHOD OF BLOCKING FOR NON-SYMMETRIC SPARSE MATRIX.         * */
/*                                                                     * */
/*  PARAMETERS:                                                        * */
/*   ON ENTRY:                                                         * */
/*     IA     THE ARRAY WHICH CONTAINS ROW INDEX OF NON-ZERO ELEMENTS  * */
/*            OF THE MATRIX.                                           * */
/*     L      THE LEADING DIMENSION OF THE ARRAY A.                    * */
/*     MA     ACCUMULATED SUM OF NUMBERS OF NON-ZERO ELEMENTS IN       * */
/*            EACH COLUMN OF THE MATRIX A.                             * */
/*     M      THE ORDER OF THE MATRIX A.                               * */
/*   ON RETURN:                                                        * */
/*     JROW   THE INFORMATION ABOUT CHANGE OF ROWS.                    * */
/*     JCOL   THE INFORMATION ABOUT CHANGE OF COLUMNS.                 * */
/*     IW     BLOCK INDEX OF EACH COLUMN.                              * */
/*     LG     NUMBER OF BLOCKS OF THE KERNEL.                          * */
/*     IER    THE ERROR CODE. IF IER=0, NORMAL RETURN.                 * */
/*   OTHERS:  WORKING PARAMETERS.                                      * */
/*                                                                     * */
/*  COPYRIGHT:     TSUTOMU OGUNI     SEP. 1 1991        VER. 1         * */
/* ********************************************************************** */

    /* Parameter adjustments */
    --ia;
    --iw;
    --jp;
    --ip;
    --jcol;
    --jrow;
    --ic;
    --ir;
    --c__;
    --r__;

    /* Function Body */
    *ier = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__[i__] = 0;
	c__[i__] = 0;
	iw[i__] = 0;
	jrow[i__] = 0;
	jcol[i__] = 0;
	ir[i__] = 0;
/* L1: */
	ic[i__] = 0;
    }

    single_(&ia[1], l, &ir[1], &ic[1], &r__[1], &c__[1], ma, m, mend, kerns, &
	    jcol[1], &jrow[1]);

    printf("m = %ld, mend = %ld, kerns = %ld\n",*m,*mend,*kerns);
    for (i__ = 1; i__ <= *l; ++i__) printf("ia[%ld]=%ld\n",i__,ia[i__]);
    for (i__ = 1; i__ <= i__1; ++i__) printf("ir[%ld]=%ld\n",i__,ir[i__]);
    for (i__ = 1; i__ <= i__1; ++i__) printf("ic[%ld]=%ld\n",i__,ic[i__]);
    for (i__ = 1; i__ <= i__1; ++i__) printf("r__[%ld]=%ld\n",i__,r__[i__]);
    for (i__ = 1; i__ <= i__1; ++i__) printf("c__[%ld]=%ld\n",i__,c__[i__]);
    for (i__ = 1; i__ <= i__1; ++i__) printf("ma[%ld]=%ld\n",i__,ma[i__]);
    for (i__ = 1; i__ <= i__1; ++i__) printf("jcol[%ld]=%ld\n",i__,jcol[i__]);
    for (i__ = 1; i__ <= i__1; ++i__) printf("jrow[%ld]=%ld\n",i__,jrow[i__]);

    forful_(&ia[1], l, ma, m, &ip[1], &jp[1], &ir[1], &ic[1], kerns, ier);

    kern = *kerns;
    mm = *m;
    *lg = 1;
L5:

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__[i__] = 0;
/* L10: */
	c__[i__] = 0;
    }
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	if (jcol[j] == 0) {
	    i__2 = ma[j];
	    for (n = ma[j - 1] + 1; n <= i__2; ++n) {
		if (jrow[ia[n]] == 0) {
		    ++r__[ia[n]];
		    ++c__[j];
		}
/* L30: */
	    }
	}
/* L20: */
    }
/*  SEARCH OF MINIMUM ROW COUNT */
    min__ = 100000;
    irow = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (jrow[i__] == 0) {
	    if (r__[i__] < min__) {
		min__ = r__[i__];
		irow = i__;
	    }
	}
/* L40: */
    }
    if (irow == 0) {
	s_wsle(&io___8);
	do_lio(&c__9, &c__1, "(SUBR. STWART) STOP AT. ", (ftnlen)24);
	do_lio(&c__3, &c__1, (char *)&kern, (ftnlen)sizeof(long));
	e_wsle();
	*ier = 1;
	return 0;
    }
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	if (jcol[j] == 0) {
	    i__2 = ma[j];
	    for (n = ma[j - 1] + 1; n <= i__2; ++n) {
		if (ia[n] == irow) {
		    iw[j] = *lg;
		    ++kern;
		}
/* L55: */
	    }
	}
/* L50: */
    }

L57:
    is = 0;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	if (iw[j] == *lg) {
	    irow = jp[j];
	    i__2 = *m;
	    for (k = 1; k <= i__2; ++k) {
		if (jcol[k] == 0) {
		    if (iw[k] == 0) {
			i__3 = ma[k];
			for (n = ma[k - 1] + 1; n <= i__3; ++n) {
			    if (ia[n] == irow) {
				iw[k] = *lg;
				++kern;
				++is;
			    }
/* L80: */
			}
		    }
		}
/* L70: */
	    }
	}
/* L60: */
    }

    if (is != 0) {
	goto L57;
    }
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	if (iw[j] == *lg) {
	    jcol[j] = mm;
	    jrow[jp[j]] = mm;
	    --mm;
	}
/* L90: */
    }
    if (*lg < *m) {
	if (kern <= *m) {
	    ++(*lg);
	    goto L5;
	}
    }
    *ier = 0;
    return 0;
/*  END OF STWART */
} /* stwart_ */

