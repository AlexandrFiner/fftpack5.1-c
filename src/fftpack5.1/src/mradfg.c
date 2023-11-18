/* f2c.h  --  Standard Fortran to C header file */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

	- From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>
#ifdef complex
#undef complex
#endif
#ifdef I
#undef I
#endif

typedef int integer;
typedef unsigned int uinteger;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
static inline _Complex float Cf(complex *z) {return z->r + z->i*_Complex_I;}
static inline _Complex double Cd(doublecomplex *z) {return z->r + z->i*_Complex_I;}
static inline _Complex float * _pCf(complex *z) {return (_Complex float*)z;}
static inline _Complex double * _pCd(doublecomplex *z) {return (_Complex double*)z;}
#define pCf(z) (*_pCf(z))
#define pCd(z) (*_pCd(z))
typedef int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

/* I/O stuff */

typedef int flag;
typedef int ftnlen;
typedef int ftnint;

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
	integer1 g;
	shortint h;
	integer i;
	/* longint j; */
	real r;
	doublereal d;
	complex c;
	doublecomplex z;
	};

typedef union Multitype Multitype;

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
#define dabs(x) (fabs(x))
#define f2cmin(a,b) ((a) <= (b) ? (a) : (b))
#define f2cmax(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (f2cmin(a,b))
#define dmax(a,b) (f2cmax(a,b))
#define bit_test(a,b)	((a) >> (b) & 1)
#define bit_clear(a,b)	((a) & ~((uinteger)1 << (b)))
#define bit_set(a,b)	((a) |  ((uinteger)1 << (b)))

#define abort_() { sig_die("Fortran abort routine called", 1); }
#define c_abs(z) (cabsf(Cf(z)))
#define c_cos(R,Z) { pCf(R)=ccos(Cf(Z)); }
#define c_div(c, a, b) {pCf(c) = Cf(a)/Cf(b);}
#define z_div(c, a, b) {pCd(c) = Cd(a)/Cd(b);}
#define c_exp(R, Z) {pCf(R) = cexpf(Cf(Z));}
#define c_log(R, Z) {pCf(R) = clogf(Cf(Z));}
#define c_sin(R, Z) {pCf(R) = csinf(Cf(Z));}
#define c_sqrt(R, Z) {*(R) = csqrtf(Cf(Z));}
#define d_abs(x) (fabs(*(x)))
#define d_acos(x) (acos(*(x)))
#define d_asin(x) (asin(*(x)))
#define d_atan(x) (atan(*(x)))
#define d_atn2(x, y) (atan2(*(x),*(y)))
#define d_cnjg(R, Z) { pCd(R) = conj(Cd(Z)); }
#define d_cos(x) (cos(*(x)))
#define d_cosh(x) (cosh(*(x)))
#define d_dim(__a, __b) ( *(__a) > *(__b) ? *(__a) - *(__b) : 0.0 )
#define d_exp(x) (exp(*(x)))
#define d_imag(z) (cimag(Cd(z)))
#define d_int(__x) (*(__x)>0 ? floor(*(__x)) : -floor(- *(__x)))
#define d_lg10(x) ( 0.43429448190325182765 * log(*(x)) )
#define d_log(x) (log(*(x)))
#define d_mod(x, y) (fmod(*(x), *(y)))
#define u_nint(__x) ((__x)>=0 ? floor((__x) + .5) : -floor(.5 - (__x)))
#define d_nint(x) u_nint(*(x))
#define u_sign(__a,__b) ((__b) >= 0 ? ((__a) >= 0 ? (__a) : -(__a)) : -((__a) >= 0 ? (__a) : -(__a)))
#define d_sign(a,b) u_sign(*(a),*(b))
#define d_sin(x) (sin(*(x)))
#define d_sinh(x) (sinh(*(x)))
#define d_sqrt(x) (sqrt(*(x)))
#define d_tan(x) (tan(*(x)))
#define d_tanh(x) (tanh(*(x)))
#define i_abs(x) abs(*(x))
#define i_dnnt(x) ((integer)u_nint(*(x)))
#define i_len(s, n) (n)
#define i_nint(x) ((integer)u_nint(*(x)))
#define i_sign(a,b) ((integer)u_sign((integer)*(a),(integer)*(b)))
#define pow_ci(p, a, b) { pCf(p) = cpow_ui(Cf(a), *(b)); }
#define pow_dd(ap, bp) ( pow(*(ap), *(bp)))
#define pow_si(B,E) spow_ui(*(B),*(E))
#define pow_di(B,E) dpow_ui(*(B),*(E))
#define pow_zi(p, a, b) {pCd(p) = zpow_ui(Cd(a), *(b));}
#define pow_zz(R,A,B) {pCd(R) = cpow(Cd(A),*(B));}
#define s_cat(lpp, rpp, rnp, np, llp) { 	ftnlen i, nc, ll; char *f__rp, *lp; 	ll = (llp); lp = (lpp); 	for(i=0; i < (int)*(np); ++i) {         	nc = ll; 	        if((rnp)[i] < nc) nc = (rnp)[i]; 	        ll -= nc;         	f__rp = (rpp)[i]; 	        while(--nc >= 0) *lp++ = *(f__rp)++;         } 	while(--ll >= 0) *lp++ = ' '; }
#define s_cmp(a,b,c,d) ((integer)strncmp((a),(b),f2cmin((c),(d))))
#define s_copy(A,B,C,D) { int __i,__m; for (__i=0, __m=f2cmin((C),(D)); __i<__m && (B)[__i] != 0; ++__i) (A)[__i] = (B)[__i]; }
#define sig_die(s, kill) { exit(1); }
#define s_stop(s, n) {exit(0);}
static char junk[] = "\n@(#)LIBF77 VERSION 19990503\n";
#define z_abs(z) (cabs(Cd(z)))
#define z_exp(R, Z) {pCd(R) = cexp(Cd(Z));}
#define z_sqrt(R, Z) {pCd(R) = csqrt(Cd(Z));}
#define myexit_() break;
#define mymaxloc_(w,s,e,n) {if (sizeof(*(w)) == sizeof(double)) dmaxloc_((w),*(s),*(e),n); else dmaxloc_((w),*(s),*(e),n);}

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef logical (*L_fp)(...);
#else
typedef logical (*L_fp)();
#endif

static float spow_ui(float x, integer n) {
	float pow=1.0; unsigned long int u;
	if(n != 0) {
		if(n < 0) n = -n, x = 1/x;
		for(u = n; ; ) {
			if(u & 01) pow *= x;
			if(u >>= 1) x *= x;
			else break;
		}
	}
	return pow;
}
static double dpow_ui(double x, integer n) {
	double pow=1.0; unsigned long int u;
	if(n != 0) {
		if(n < 0) n = -n, x = 1/x;
		for(u = n; ; ) {
			if(u & 01) pow *= x;
			if(u >>= 1) x *= x;
			else break;
		}
	}
	return pow;
}
static _Complex float cpow_ui(_Complex float x, integer n) {
	_Complex float pow=1.0; unsigned long int u;
	if(n != 0) {
		if(n < 0) n = -n, x = 1/x;
		for(u = n; ; ) {
			if(u & 01) pow *= x;
			if(u >>= 1) x *= x;
			else break;
		}
	}
	return pow;
}
static _Complex double zpow_ui(_Complex double x, integer n) {
	_Complex double pow=1.0; unsigned long int u;
	if(n != 0) {
		if(n < 0) n = -n, x = 1/x;
		for(u = n; ; ) {
			if(u & 01) pow *= x;
			if(u >>= 1) x *= x;
			else break;
		}
	}
	return pow;
}
static integer pow_ii(integer x, integer n) {
	integer pow; unsigned long int u;
	if (n <= 0) {
		if (n == 0 || x == 1) pow = 1;
		else if (x != -1) pow = x == 0 ? 1/x : 0;
		else n = -n;
	}
	if ((n > 0) || !(n == 0 || x == 1 || x != -1)) {
		u = n;
		for(pow = 1; ; ) {
			if(u & 01) pow *= x;
			if(u >>= 1) x *= x;
			else break;
		}
	}
	return pow;
}
static integer dmaxloc_(double *w, integer s, integer e, integer *n)
{
	double m; integer i, mi;
	for(m=w[s-1], mi=s, i=s+1; i<=e; i++)
		if (w[i-1]>m) mi=i ,m=w[i-1];
	return mi-s+1;
}
static integer smaxloc_(float *w, integer s, integer e, integer *n)
{
	float m; integer i, mi;
	for(m=w[s-1], mi=s, i=s+1; i<=e; i++)
		if (w[i-1]>m) mi=i ,m=w[i-1];
	return mi-s+1;
}
static inline void cdotc_(complex *z, integer *n_, complex *x, integer *incx_, complex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
	_Complex float zdotc = 0.0;
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += conjf(Cf(&x[i])) * Cf(&y[i]);
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += conjf(Cf(&x[i*incx])) * Cf(&y[i*incy]);
		}
	}
	pCf(z) = zdotc;
}
static inline void zdotc_(doublecomplex *z, integer *n_, doublecomplex *x, integer *incx_, doublecomplex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
	_Complex double zdotc = 0.0;
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += conj(Cd(&x[i])) * Cd(&y[i]);
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += conj(Cd(&x[i*incx])) * Cd(&y[i*incy]);
		}
	}
	pCd(z) = zdotc;
}
static inline void cdotu_(complex *z, integer *n_, complex *x, integer *incx_, complex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
	_Complex float zdotc = 0.0;
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += Cf(&x[i]) * Cf(&y[i]);
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += Cf(&x[i*incx]) * Cf(&y[i*incy]);
		}
	}
	pCf(z) = zdotc;
}
static inline void zdotu_(doublecomplex *z, integer *n_, doublecomplex *x, integer *incx_, doublecomplex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
	_Complex double zdotc = 0.0;
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += Cd(&x[i]) * Cd(&y[i]);
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += Cd(&x[i*incx]) * Cd(&y[i*incy]);
		}
	}
	pCd(z) = zdotc;
}
#endif
/*  -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*     *                                                               * */
/*     *                  copyright (c) 2011 by UCAR                   * */
/*     *                                                               * */
/*     *       University Corporation for Atmospheric Research         * */
/*     *                                                               * */
/*     *                      all rights reserved                      * */
/*     *                                                               * */
/*     *                     FFTPACK  version 5.1                      * */
/*     *                                                               * */
/*     *                 A Fortran Package of Fast Fourier             * */
/*     *                                                               * */
/*     *                Subroutines and Example Programs               * */
/*     *                                                               * */
/*     *                             by                                * */
/*     *                                                               * */
/*     *               Paul Swarztrauber and Dick Valent               * */
/*     *                                                               * */
/*     *                             of                                * */
/*     *                                                               * */
/*     *         the National Center for Atmospheric Research          * */
/*     *                                                               * */
/*     *                Boulder, Colorado  (80307)  U.S.A.             * */
/*     *                                                               * */
/*     *                   which is sponsored by                       * */
/*     *                                                               * */
/*     *              the National Science Foundation                  * */
/*     *                                                               * */
/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Subroutine */ int mradfg_(integer *m, integer *ido, integer *ip, integer *
	l1, integer *idl1, real *cc, real *c1, real *c2, integer *im1, 
	integer *in1, real *ch, real *ch2, integer *im2, integer *in2, real *
	wa)
{
    /* System generated locals */
    integer ch_dim1, ch_dim2, ch_dim3, ch_offset, cc_dim1, cc_dim2, cc_dim3, 
	    cc_offset, c1_dim1, c1_dim2, c1_dim3, c1_offset, c2_dim1, c2_dim2,
	     c2_offset, ch2_dim1, ch2_dim2, ch2_offset, i__1, i__2, i__3, 
	    i__4, i__5;

    /* Local variables */
    integer i__, j, k, l, j2, m1, m2, ic, jc, lc, ik, is;
    real dc2, ai1, ai2;
    integer m1d;
    real ar1, ar2, ds2;
    integer m2s, nbd;
    real dcp, arg, dsp, tpi, ar1h, ar2h;
    integer idp2, ipp2, idij, ipph;


    /* Parameter adjustments */
    --wa;
    c2_dim1 = *in1;
    c2_dim2 = *idl1;
    c2_offset = 1 + c2_dim1 * (1 + c2_dim2);
    c2 -= c2_offset;
    c1_dim1 = *in1;
    c1_dim2 = *ido;
    c1_dim3 = *l1;
    c1_offset = 1 + c1_dim1 * (1 + c1_dim2 * (1 + c1_dim3));
    c1 -= c1_offset;
    cc_dim1 = *in1;
    cc_dim2 = *ido;
    cc_dim3 = *ip;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * (1 + cc_dim3));
    cc -= cc_offset;
    ch2_dim1 = *in2;
    ch2_dim2 = *idl1;
    ch2_offset = 1 + ch2_dim1 * (1 + ch2_dim2);
    ch2 -= ch2_offset;
    ch_dim1 = *in2;
    ch_dim2 = *ido;
    ch_dim3 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * (1 + ch_dim3));
    ch -= ch_offset;

    /* Function Body */
    m1d = (*m - 1) * *im1 + 1;
    m2s = 1 - *im2;
    tpi = atan(1.f) * 8.f;
    arg = tpi / (real) (*ip);
    dcp = cos(arg);
    dsp = sin(arg);
    ipph = (*ip + 1) / 2;
    ipp2 = *ip + 2;
    idp2 = *ido + 2;
    nbd = (*ido - 1) / 2;
    if (*ido == 1) {
	goto L119;
    }
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	m2 = m2s;
	i__2 = m1d;
	i__3 = *im1;
	for (m1 = 1; i__3 < 0 ? m1 >= i__2 : m1 <= i__2; m1 += i__3) {
	    m2 += *im2;
	    ch2[m2 + (ik + ch2_dim2) * ch2_dim1] = c2[m1 + (ik + c2_dim2) * 
		    c2_dim1];
/* L1001: */
	}
/* L101: */
    }
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	i__3 = *l1;
	for (k = 1; k <= i__3; ++k) {
	    m2 = m2s;
	    i__2 = m1d;
	    i__4 = *im1;
	    for (m1 = 1; i__4 < 0 ? m1 >= i__2 : m1 <= i__2; m1 += i__4) {
		m2 += *im2;
		ch[m2 + ((k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1] = c1[m1 
			+ ((k + j * c1_dim3) * c1_dim2 + 1) * c1_dim1];
/* L1002: */
	    }
/* L102: */
	}
/* L103: */
    }
    if (nbd > *l1) {
	goto L107;
    }
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	idij = is;
	i__3 = *ido;
	for (i__ = 3; i__ <= i__3; i__ += 2) {
	    idij += 2;
	    i__4 = *l1;
	    for (k = 1; k <= i__4; ++k) {
		m2 = m2s;
		i__2 = m1d;
		i__5 = *im1;
		for (m1 = 1; i__5 < 0 ? m1 >= i__2 : m1 <= i__2; m1 += i__5) {
		    m2 += *im2;
		    ch[m2 + (i__ - 1 + (k + j * ch_dim3) * ch_dim2) * ch_dim1]
			     = wa[idij - 1] * c1[m1 + (i__ - 1 + (k + j * 
			    c1_dim3) * c1_dim2) * c1_dim1] + wa[idij] * c1[m1 
			    + (i__ + (k + j * c1_dim3) * c1_dim2) * c1_dim1];
		    ch[m2 + (i__ + (k + j * ch_dim3) * ch_dim2) * ch_dim1] = 
			    wa[idij - 1] * c1[m1 + (i__ + (k + j * c1_dim3) * 
			    c1_dim2) * c1_dim1] - wa[idij] * c1[m1 + (i__ - 1 
			    + (k + j * c1_dim3) * c1_dim2) * c1_dim1];
/* L1004: */
		}
/* L104: */
	    }
/* L105: */
	}
/* L106: */
    }
    goto L111;
L107:
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	i__3 = *l1;
	for (k = 1; k <= i__3; ++k) {
	    idij = is;
	    i__4 = *ido;
	    for (i__ = 3; i__ <= i__4; i__ += 2) {
		idij += 2;
		m2 = m2s;
		i__5 = m1d;
		i__2 = *im1;
		for (m1 = 1; i__2 < 0 ? m1 >= i__5 : m1 <= i__5; m1 += i__2) {
		    m2 += *im2;
		    ch[m2 + (i__ - 1 + (k + j * ch_dim3) * ch_dim2) * ch_dim1]
			     = wa[idij - 1] * c1[m1 + (i__ - 1 + (k + j * 
			    c1_dim3) * c1_dim2) * c1_dim1] + wa[idij] * c1[m1 
			    + (i__ + (k + j * c1_dim3) * c1_dim2) * c1_dim1];
		    ch[m2 + (i__ + (k + j * ch_dim3) * ch_dim2) * ch_dim1] = 
			    wa[idij - 1] * c1[m1 + (i__ + (k + j * c1_dim3) * 
			    c1_dim2) * c1_dim1] - wa[idij] * c1[m1 + (i__ - 1 
			    + (k + j * c1_dim3) * c1_dim2) * c1_dim1];
/* L1008: */
		}
/* L108: */
	    }
/* L109: */
	}
/* L110: */
    }
L111:
    if (nbd < *l1) {
	goto L115;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__3 = *l1;
	for (k = 1; k <= i__3; ++k) {
	    i__4 = *ido;
	    for (i__ = 3; i__ <= i__4; i__ += 2) {
		m2 = m2s;
		i__2 = m1d;
		i__5 = *im1;
		for (m1 = 1; i__5 < 0 ? m1 >= i__2 : m1 <= i__2; m1 += i__5) {
		    m2 += *im2;
		    c1[m1 + (i__ - 1 + (k + j * c1_dim3) * c1_dim2) * c1_dim1]
			     = ch[m2 + (i__ - 1 + (k + j * ch_dim3) * ch_dim2)
			     * ch_dim1] + ch[m2 + (i__ - 1 + (k + jc * 
			    ch_dim3) * ch_dim2) * ch_dim1];
		    c1[m1 + (i__ - 1 + (k + jc * c1_dim3) * c1_dim2) * 
			    c1_dim1] = ch[m2 + (i__ + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] - ch[m2 + (i__ + (k + jc * 
			    ch_dim3) * ch_dim2) * ch_dim1];
		    c1[m1 + (i__ + (k + j * c1_dim3) * c1_dim2) * c1_dim1] = 
			    ch[m2 + (i__ + (k + j * ch_dim3) * ch_dim2) * 
			    ch_dim1] + ch[m2 + (i__ + (k + jc * ch_dim3) * 
			    ch_dim2) * ch_dim1];
		    c1[m1 + (i__ + (k + jc * c1_dim3) * c1_dim2) * c1_dim1] = 
			    ch[m2 + (i__ - 1 + (k + jc * ch_dim3) * ch_dim2) *
			     ch_dim1] - ch[m2 + (i__ - 1 + (k + j * ch_dim3) *
			     ch_dim2) * ch_dim1];
/* L1012: */
		}
/* L112: */
	    }
/* L113: */
	}
/* L114: */
    }
    goto L121;
L115:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__3 = *ido;
	for (i__ = 3; i__ <= i__3; i__ += 2) {
	    i__4 = *l1;
	    for (k = 1; k <= i__4; ++k) {
		m2 = m2s;
		i__5 = m1d;
		i__2 = *im1;
		for (m1 = 1; i__2 < 0 ? m1 >= i__5 : m1 <= i__5; m1 += i__2) {
		    m2 += *im2;
		    c1[m1 + (i__ - 1 + (k + j * c1_dim3) * c1_dim2) * c1_dim1]
			     = ch[m2 + (i__ - 1 + (k + j * ch_dim3) * ch_dim2)
			     * ch_dim1] + ch[m2 + (i__ - 1 + (k + jc * 
			    ch_dim3) * ch_dim2) * ch_dim1];
		    c1[m1 + (i__ - 1 + (k + jc * c1_dim3) * c1_dim2) * 
			    c1_dim1] = ch[m2 + (i__ + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] - ch[m2 + (i__ + (k + jc * 
			    ch_dim3) * ch_dim2) * ch_dim1];
		    c1[m1 + (i__ + (k + j * c1_dim3) * c1_dim2) * c1_dim1] = 
			    ch[m2 + (i__ + (k + j * ch_dim3) * ch_dim2) * 
			    ch_dim1] + ch[m2 + (i__ + (k + jc * ch_dim3) * 
			    ch_dim2) * ch_dim1];
		    c1[m1 + (i__ + (k + jc * c1_dim3) * c1_dim2) * c1_dim1] = 
			    ch[m2 + (i__ - 1 + (k + jc * ch_dim3) * ch_dim2) *
			     ch_dim1] - ch[m2 + (i__ - 1 + (k + j * ch_dim3) *
			     ch_dim2) * ch_dim1];
/* L1016: */
		}
/* L116: */
	    }
/* L117: */
	}
/* L118: */
    }
    goto L121;
L119:
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	m2 = m2s;
	i__3 = m1d;
	i__4 = *im1;
	for (m1 = 1; i__4 < 0 ? m1 >= i__3 : m1 <= i__3; m1 += i__4) {
	    m2 += *im2;
	    c2[m1 + (ik + c2_dim2) * c2_dim1] = ch2[m2 + (ik + ch2_dim2) * 
		    ch2_dim1];
/* L1020: */
	}
/* L120: */
    }
L121:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__4 = *l1;
	for (k = 1; k <= i__4; ++k) {
	    m2 = m2s;
	    i__3 = m1d;
	    i__2 = *im1;
	    for (m1 = 1; i__2 < 0 ? m1 >= i__3 : m1 <= i__3; m1 += i__2) {
		m2 += *im2;
		c1[m1 + ((k + j * c1_dim3) * c1_dim2 + 1) * c1_dim1] = ch[m2 
			+ ((k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1] + ch[
			m2 + ((k + jc * ch_dim3) * ch_dim2 + 1) * ch_dim1];
		c1[m1 + ((k + jc * c1_dim3) * c1_dim2 + 1) * c1_dim1] = ch[m2 
			+ ((k + jc * ch_dim3) * ch_dim2 + 1) * ch_dim1] - ch[
			m2 + ((k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1];
/* L1022: */
	    }
/* L122: */
	}
/* L123: */
    }

    ar1 = 1.f;
    ai1 = 0.f;
    i__1 = ipph;
    for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	ar1h = dcp * ar1 - dsp * ai1;
	ai1 = dcp * ai1 + dsp * ar1;
	ar1 = ar1h;
	i__4 = *idl1;
	for (ik = 1; ik <= i__4; ++ik) {
	    m2 = m2s;
	    i__2 = m1d;
	    i__3 = *im1;
	    for (m1 = 1; i__3 < 0 ? m1 >= i__2 : m1 <= i__2; m1 += i__3) {
		m2 += *im2;
		ch2[m2 + (ik + l * ch2_dim2) * ch2_dim1] = c2[m1 + (ik + 
			c2_dim2) * c2_dim1] + ar1 * c2[m1 + (ik + (c2_dim2 << 
			1)) * c2_dim1];
		ch2[m2 + (ik + lc * ch2_dim2) * ch2_dim1] = ai1 * c2[m1 + (ik 
			+ *ip * c2_dim2) * c2_dim1];
/* L1024: */
	    }
/* L124: */
	}
	dc2 = ar1;
	ds2 = ai1;
	ar2 = ar1;
	ai2 = ai1;
	i__4 = ipph;
	for (j = 3; j <= i__4; ++j) {
	    jc = ipp2 - j;
	    ar2h = dc2 * ar2 - ds2 * ai2;
	    ai2 = dc2 * ai2 + ds2 * ar2;
	    ar2 = ar2h;
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		m2 = m2s;
		i__2 = m1d;
		i__5 = *im1;
		for (m1 = 1; i__5 < 0 ? m1 >= i__2 : m1 <= i__2; m1 += i__5) {
		    m2 += *im2;
		    ch2[m2 + (ik + l * ch2_dim2) * ch2_dim1] += ar2 * c2[m1 + 
			    (ik + j * c2_dim2) * c2_dim1];
		    ch2[m2 + (ik + lc * ch2_dim2) * ch2_dim1] += ai2 * c2[m1 
			    + (ik + jc * c2_dim2) * c2_dim1];
/* L1025: */
		}
/* L125: */
	    }
/* L126: */
	}
/* L127: */
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	i__4 = *idl1;
	for (ik = 1; ik <= i__4; ++ik) {
	    m2 = m2s;
	    i__3 = m1d;
	    i__5 = *im1;
	    for (m1 = 1; i__5 < 0 ? m1 >= i__3 : m1 <= i__3; m1 += i__5) {
		m2 += *im2;
		ch2[m2 + (ik + ch2_dim2) * ch2_dim1] += c2[m1 + (ik + j * 
			c2_dim2) * c2_dim1];
/* L1028: */
	    }
/* L128: */
	}
/* L129: */
    }

    if (*ido < *l1) {
	goto L132;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__4 = *ido;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    m2 = m2s;
	    i__5 = m1d;
	    i__3 = *im1;
	    for (m1 = 1; i__3 < 0 ? m1 >= i__5 : m1 <= i__5; m1 += i__3) {
		m2 += *im2;
		cc[m1 + (i__ + (k * cc_dim3 + 1) * cc_dim2) * cc_dim1] = ch[
			m2 + (i__ + (k + ch_dim3) * ch_dim2) * ch_dim1];
/* L1030: */
	    }
/* L130: */
	}
/* L131: */
    }
    goto L135;
L132:
    i__1 = *ido;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__4 = *l1;
	for (k = 1; k <= i__4; ++k) {
	    m2 = m2s;
	    i__3 = m1d;
	    i__5 = *im1;
	    for (m1 = 1; i__5 < 0 ? m1 >= i__3 : m1 <= i__3; m1 += i__5) {
		m2 += *im2;
		cc[m1 + (i__ + (k * cc_dim3 + 1) * cc_dim2) * cc_dim1] = ch[
			m2 + (i__ + (k + ch_dim3) * ch_dim2) * ch_dim1];
/* L1033: */
	    }
/* L133: */
	}
/* L134: */
    }
L135:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__4 = *l1;
	for (k = 1; k <= i__4; ++k) {
	    m2 = m2s;
	    i__5 = m1d;
	    i__3 = *im1;
	    for (m1 = 1; i__3 < 0 ? m1 >= i__5 : m1 <= i__5; m1 += i__3) {
		m2 += *im2;
		cc[m1 + (*ido + (j2 - 2 + k * cc_dim3) * cc_dim2) * cc_dim1] =
			 ch[m2 + ((k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1];
		cc[m1 + ((j2 - 1 + k * cc_dim3) * cc_dim2 + 1) * cc_dim1] = 
			ch[m2 + ((k + jc * ch_dim3) * ch_dim2 + 1) * ch_dim1];
/* L1036: */
	    }
/* L136: */
	}
/* L137: */
    }
    if (*ido == 1) {
	return 0;
    }
    if (nbd < *l1) {
	goto L141;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__4 = *l1;
	for (k = 1; k <= i__4; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		ic = idp2 - i__;
		m2 = m2s;
		i__5 = m1d;
		i__2 = *im1;
		for (m1 = 1; i__2 < 0 ? m1 >= i__5 : m1 <= i__5; m1 += i__2) {
		    m2 += *im2;
		    cc[m1 + (i__ - 1 + (j2 - 1 + k * cc_dim3) * cc_dim2) * 
			    cc_dim1] = ch[m2 + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] + ch[m2 + (i__ - 1 + (k + jc *
			     ch_dim3) * ch_dim2) * ch_dim1];
		    cc[m1 + (ic - 1 + (j2 - 2 + k * cc_dim3) * cc_dim2) * 
			    cc_dim1] = ch[m2 + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] - ch[m2 + (i__ - 1 + (k + jc *
			     ch_dim3) * ch_dim2) * ch_dim1];
		    cc[m1 + (i__ + (j2 - 1 + k * cc_dim3) * cc_dim2) * 
			    cc_dim1] = ch[m2 + (i__ + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] + ch[m2 + (i__ + (k + jc * 
			    ch_dim3) * ch_dim2) * ch_dim1];
		    cc[m1 + (ic + (j2 - 2 + k * cc_dim3) * cc_dim2) * cc_dim1]
			     = ch[m2 + (i__ + (k + jc * ch_dim3) * ch_dim2) * 
			    ch_dim1] - ch[m2 + (i__ + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1];
/* L1038: */
		}
/* L138: */
	    }
/* L139: */
	}
/* L140: */
    }
    return 0;
L141:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__4 = *ido;
	for (i__ = 3; i__ <= i__4; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		m2 = m2s;
		i__2 = m1d;
		i__5 = *im1;
		for (m1 = 1; i__5 < 0 ? m1 >= i__2 : m1 <= i__2; m1 += i__5) {
		    m2 += *im2;
		    cc[m1 + (i__ - 1 + (j2 - 1 + k * cc_dim3) * cc_dim2) * 
			    cc_dim1] = ch[m2 + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] + ch[m2 + (i__ - 1 + (k + jc *
			     ch_dim3) * ch_dim2) * ch_dim1];
		    cc[m1 + (ic - 1 + (j2 - 2 + k * cc_dim3) * cc_dim2) * 
			    cc_dim1] = ch[m2 + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] - ch[m2 + (i__ - 1 + (k + jc *
			     ch_dim3) * ch_dim2) * ch_dim1];
		    cc[m1 + (i__ + (j2 - 1 + k * cc_dim3) * cc_dim2) * 
			    cc_dim1] = ch[m2 + (i__ + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] + ch[m2 + (i__ + (k + jc * 
			    ch_dim3) * ch_dim2) * ch_dim1];
		    cc[m1 + (ic + (j2 - 2 + k * cc_dim3) * cc_dim2) * cc_dim1]
			     = ch[m2 + (i__ + (k + jc * ch_dim3) * ch_dim2) * 
			    ch_dim1] - ch[m2 + (i__ + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1];
/* L1042: */
		}
/* L142: */
	    }
/* L143: */
	}
/* L144: */
    }
    return 0;
} /* mradfg_ */

