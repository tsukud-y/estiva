/* SINGLE.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int single_(ia, l, ir, ic, r__, c__, ma, m, mend, kerns, 
	jcol, jrow)
integer *ia, *l, *ir, *ic, *r__, *c__, *ma, *m, *mend, *kerns, *jcol, *jrow;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer icol, irow, i__, j, k, n;

/* ********************************************************************** */
/*  DETECTION OF SINGLETONS OF A SPARSE MATRIX.                        * */
/*                                                                     * */
/*  PARAMETERS:                                                        * */
/*   ON ENTRY:                                                         * */
/*     SAME AS MARKWZ ROUTINE.                                         * */
/*   ON RETURN:                                                        * */
/*     IR     ROW INDEX OF ROW SINGLETONS IN THE FIRST MEND ENTRIES    * */
/*            IN THE ARRAY IR. ROW INDEX OF COLUMN SINGLETONS IN       * */
/*            (KERNS-MEND-1) ENTRIES IN THE ARRAY IR.                  * */
/*     IC     COLUMN INDEX OF ROW SINGLETONS IN THE FIRST MEND ENTRIES * */
/*            IN THE ARRAY IC. COLUMN INDEX OF COLUMN SINGLETON IN     * */
/*            (KERNS-MEND-1) ENTRIES IN THE ARRAY IC.                  * */
/*   OTHERS:  THE WORKING PARAMETERS.                                  * */
/*                                                                     * */
/*  COPYRIGHT:      TSUTOMU OGUNI   SEP. 1 1992       VER. 2           * */
/* ********************************************************************** */
/*  SINGLETON */
/*  ROW SINGLETON */
    /* Parameter adjustments */
    --ia;
    --jrow;
    --jcol;
    --c__;
    --r__;
    --ic;
    --ir;

    /* Function Body */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	jrow[i__] = 0;
	jcol[i__] = 0;
	ic[i__] = 0;
	ir[i__] = 0;
	c__[i__] = ma[i__] - ma[i__ - 1];
	i__2 = ma[i__];
	for (n = ma[i__ - 1] + 1; n <= i__2; ++n) {
/* L110: */
	    ++r__[ia[n]];
	}
/* L100: */
    }
    *kerns = 1;
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	irow = 0;
	icol = 0;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (r__[i__] == 1) {
		irow = i__;
		goto L170;
	    }
/* L160: */
	}
L170:
	if (irow == 0) {
	    goto L180;
	}
	ir[*kerns] = irow;
	jrow[irow] = *kerns;

	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    if (jcol[j] == 0) {
		i__3 = ma[j];
		for (n = ma[j - 1] + 1; n <= i__3; ++n) {
		    if (ia[n] == irow) {
			ic[*kerns] = j;
			icol = j;
			goto L200;
		    }
/* L190: */
		}
	    }
L200:
	    ;
	}
	jcol[icol] = *kerns;
	++(*kerns);
	i__2 = ma[icol];
	for (n = ma[icol - 1] + 1; n <= i__2; ++n) {
	    --r__[ia[n]];
/* L210: */
	}
/* L150: */
    }

L180:
    *mend = *kerns - 1;
/*  COLUMN SINGLETON */
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	irow = 0;
	icol = 0;
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    if (jcol[j] == 0) {
		if (c__[j] == 1) {
		    icol = j;
		    goto L135;
		}
	    }
/* L130: */
	}
L135:
	if (icol == 0) {
	    return 0;
	}
	jcol[icol] = *kerns;
	ic[*kerns] = icol;

	i__2 = ma[icol];
	for (n = ma[icol - 1] + 1; n <= i__2; ++n) {
	    if (jrow[ia[n]] == 0) {
		irow = ia[n];
	    }
/* L195: */
	}
	ir[*kerns] = irow;
	jrow[irow] = *kerns;
	++(*kerns);
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    if (jcol[j] == 0) {
		i__3 = ma[j];
		for (n = ma[j - 1] + 1; n <= i__3; ++n) {
		    if (ia[n] == irow) {
			--c__[j];
		    }
/* L185: */
		}
	    }
/* L175: */
	}
/* L120: */
    }

    return 0;
/*  END OF SINGLE */
} /* single_ */

