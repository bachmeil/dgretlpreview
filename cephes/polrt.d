/* Port of cephes/polrt.c: the complex roots of a real polynomial of degree m (at most 36), polrt(xcof, cof, m,
   root). xcof holds the m + 1 coefficients in ascending order, cof is working storage of the same size, and the
   roots are returned in root (complex roots in conjugate pairs). Newton's method on a polynomial that is deflated
   as roots are found, with a final pass on the original polynomial to polish each root. Returns 0, or 1 if m <= 0,
   2 if m > 36, 3 if the iteration fails and 4 if the leading coefficient is zero.

   Cephes Math Library, Copyright 1984, 1987, 1989 by Stephen L. Moshier.

   The goto structure of the C (nxtrut, tryagn, finitr, lupdon, findon, zerrut, setrut) and its pointer walking are
   kept, including the jump into the else block at zerrut. The variable named final in the C is final_. The
   complex type is cephes.mconfh.Cmplx (the cmplx of mconf.h); it has the layout of libgretlh.Cmplx, which
   gretl_matrix.d passes through a cast pointer. Only gretl_matrix.d calls polrt, so it has D linkage. */

module cephes.polrt;

import core.stdc.math : fabs;

import cephes.mconfh : Cmplx;

int polrt(double * xcof, double * cof, int m, Cmplx * root)
{
    double * p;
    double * q;
    int i, j, nsav, n, n1, n2, nroot, iter, retry;
    int final_;
    double mag, cofj;
    Cmplx x0, x, xsav, dx, t, t1, u, ud;

    final_ = 0;
    n = m;

    if (n <= 0) {
	return 1;
    }
    if (n > 36) {
	return 2;
    }
    if (xcof[m] == 0.0) {
	return 4;
    }

    n1 = n;
    n2 = n;
    nroot = 0;
    nsav = n;
    q = &xcof[0];
    p = &cof[n];
    for (j=0; j<=nsav; j++) {
	*p-- = *q++; /*	cof[n-j] = xcof[j]; */
    }
    xsav.r = 0.0;
    xsav.i = 0.0;

 nxtrut:
    x0.r = 0.00500101;
    x0.i = 0.01000101;
    retry = 0;

 tryagn:
    retry += 1;
    x.r = x0.r;

    x0.r = -10.0 * x0.i;
    x0.i = -10.0 * x.r;

    x.r = x0.r;
    x.i = x0.i;

 finitr:
    iter = 0;

    while (iter < 500) {
	u.r = cof[n];
	if (u.r == 0.0) {
	    /* this root is zero */
	    x.r = 0;
	    n1 -= 1;
	    n2 -= 1;
	    goto zerrut;
	}
	u.i = 0;
	ud.r = 0;
	ud.i = 0;
	t.r = 1.0;
	t.i = 0;
	p = &cof[n-1];
	for (i=0; i<n; i++) {
	    t1.r = x.r * t.r - x.i * t.i;
	    t1.i = x.r * t.i + x.i * t.r;
	    cofj = *p--;         /* evaluate polynomial */
	    u.r += cofj * t1.r;
	    u.i += cofj * t1.i;
	    cofj = cofj * (i+1); /* derivative */
	    ud.r += cofj * t.r;
	    ud.i -= cofj * t.i;
	    t.r = t1.r;
	    t.i = t1.i;
	}

	mag = ud.r * ud.r + ud.i * ud.i;
	if (mag == 0.0) {
	    if (!final_) {
		if (retry < 50) {
		    /* questionable? (this goto was not
		       conditional on the retry value, 
		       changed 2014-04-11)
		    */
		    goto tryagn;
		} else {
		    break;
		}
	    }
	    x.r = xsav.r;
	    x.i = xsav.i;
	    goto findon;
	}
	dx.r = (u.i * ud.i  -  u.r * ud.r)/mag;
	x.r += dx.r;
	dx.i = -(u.r * ud.i  +  u.i * ud.r)/mag;
	x.i += dx.i;
	if ((fabs(dx.i) + fabs(dx.r)) < 1.0e-6)
	    goto lupdon;
	iter += 1;
    } /* while iter < 500 */

    if (final_)
	goto lupdon;
    if (retry < 5)
	goto tryagn;
    return 3;

 lupdon:
    /* Swap original and reduced polynomials */
    q = &xcof[nsav];
    p = &cof[0];
    for (j=0; j<=n2; j++) {
	cofj = *q;
	*q-- = *p;
	*p++ = cofj;
    }
    i = n;
    n = n1;
    n1 = i;

    if (!final_) {
	final_ = 1;
	if(fabs(x.i/x.r) < 1.0e-4)
	    x.i = 0.0;
	xsav.r = x.r;
	xsav.i = x.i;
	goto finitr; /* do final_ iteration on original polynomial */
    }

 findon:
    final_ = 0;
    if (fabs(x.i/x.r) >= 1.0e-5) {
	cofj = x.r + x.r;
	mag = x.r * x.r  +  x.i * x.i;
	n -= 2;
    } else {
	/* root is real */
    zerrut:
	x.i = 0;
	cofj = x.r;
	mag = 0;
	n -= 1;
    }

    /* divide working polynomial cof(z) by z - x */
    p = &cof[1];
    *p += cofj * *(p-1);
    for (j=1; j<n; j++) {
	*(p+1) += cofj * *p  -  mag * *(p-1);
	p++;
    }

 setrut:

    root[nroot].r = x.r;
    root[nroot].i = x.i;
    nroot += 1;

    if (mag != 0.0) {
	x.i = -x.i;
	mag = 0;
	goto setrut;	/* fill in the complex conjugate root */
    }

    if (n > 0)
	goto nxtrut;

    return 0;
}
