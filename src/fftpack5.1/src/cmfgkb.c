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

/* Subroutine */ int cmfgkb_(integer *lot, integer *ido, integer *ip, integer 
	*l1, integer *lid, integer *na, real *cc, real *cc1, integer *im1, 
	integer *in1, real *ch, real *ch1, integer *im2, integer *in2, real *
	wa)
{
    /* System generated locals */
    integer ch_dim2, ch_dim3, ch_dim4, ch_offset, cc_dim2, cc_dim3, cc_dim4, 
	    cc_offset, cc1_dim2, cc1_dim3, cc1_offset, ch1_dim2, ch1_dim3, 
	    ch1_offset, wa_dim1, wa_dim2, wa_offset, i__1, i__2, i__3, i__4, 
	    i__5;

    /* Local variables */
    integer i__, j, k, l, m1, m2, jc, lc, ki, m1d, m2s;
    real wai, war;
    integer ipp2, idlj, ipph;
    real chold1, chold2;


/* FFTPACK 5.0 auxiliary routine */

    /* Parameter adjustments */
    wa_dim1 = *ido;
    wa_dim2 = *ip - 1;
    wa_offset = 1 + wa_dim1 * (1 + wa_dim2);
    wa -= wa_offset;
    cc1_dim2 = *in1;
    cc1_dim3 = *lid;
    cc1_offset = 1 + 2 * (1 + cc1_dim2 * (1 + cc1_dim3));
    cc1 -= cc1_offset;
    cc_dim2 = *in1;
    cc_dim3 = *l1;
    cc_dim4 = *ip;
    cc_offset = 1 + 2 * (1 + cc_dim2 * (1 + cc_dim3 * (1 + cc_dim4)));
    cc -= cc_offset;
    ch1_dim2 = *in2;
    ch1_dim3 = *lid;
    ch1_offset = 1 + 2 * (1 + ch1_dim2 * (1 + ch1_dim3));
    ch1 -= ch1_offset;
    ch_dim2 = *in2;
    ch_dim3 = *l1;
    ch_dim4 = *ido;
    ch_offset = 1 + 2 * (1 + ch_dim2 * (1 + ch_dim3 * (1 + ch_dim4)));
    ch -= ch_offset;

    /* Function Body */
    m1d = (*lot - 1) * *im1 + 1;
    m2s = 1 - *im2;
    ipp2 = *ip + 2;
    ipph = (*ip + 1) / 2;
    i__1 = *lid;
    for (ki = 1; ki <= i__1; ++ki) {
	m2 = m2s;
	i__2 = m1d;
	i__3 = *im1;
	for (m1 = 1; i__3 < 0 ? m1 >= i__2 : m1 <= i__2; m1 += i__3) {
	    m2 += *im2;
	    ch1[(m2 + (ki + ch1_dim3) * ch1_dim2 << 1) + 1] = cc1[(m1 + (ki + 
		    cc1_dim3) * cc1_dim2 << 1) + 1];
	    ch1[(m2 + (ki + ch1_dim3) * ch1_dim2 << 1) + 2] = cc1[(m1 + (ki + 
		    cc1_dim3) * cc1_dim2 << 1) + 2];
/* L110: */
	}
    }
    i__3 = ipph;
    for (j = 2; j <= i__3; ++j) {
	jc = ipp2 - j;
	i__2 = *lid;
	for (ki = 1; ki <= i__2; ++ki) {
	    m2 = m2s;
	    i__1 = m1d;
	    i__4 = *im1;
	    for (m1 = 1; i__4 < 0 ? m1 >= i__1 : m1 <= i__1; m1 += i__4) {
		m2 += *im2;
		ch1[(m2 + (ki + j * ch1_dim3) * ch1_dim2 << 1) + 1] = cc1[(m1 
			+ (ki + j * cc1_dim3) * cc1_dim2 << 1) + 1] + cc1[(m1 
			+ (ki + jc * cc1_dim3) * cc1_dim2 << 1) + 1];
		ch1[(m2 + (ki + jc * ch1_dim3) * ch1_dim2 << 1) + 1] = cc1[(
			m1 + (ki + j * cc1_dim3) * cc1_dim2 << 1) + 1] - cc1[(
			m1 + (ki + jc * cc1_dim3) * cc1_dim2 << 1) + 1];
		ch1[(m2 + (ki + j * ch1_dim3) * ch1_dim2 << 1) + 2] = cc1[(m1 
			+ (ki + j * cc1_dim3) * cc1_dim2 << 1) + 2] + cc1[(m1 
			+ (ki + jc * cc1_dim3) * cc1_dim2 << 1) + 2];
		ch1[(m2 + (ki + jc * ch1_dim3) * ch1_dim2 << 1) + 2] = cc1[(
			m1 + (ki + j * cc1_dim3) * cc1_dim2 << 1) + 2] - cc1[(
			m1 + (ki + jc * cc1_dim3) * cc1_dim2 << 1) + 2];
/* L112: */
	    }
	}
/* L111: */
    }
    i__3 = ipph;
    for (j = 2; j <= i__3; ++j) {
	i__4 = *lid;
	for (ki = 1; ki <= i__4; ++ki) {
	    m2 = m2s;
	    i__1 = m1d;
	    i__2 = *im1;
	    for (m1 = 1; i__2 < 0 ? m1 >= i__1 : m1 <= i__1; m1 += i__2) {
		m2 += *im2;
		cc1[(m1 + (ki + cc1_dim3) * cc1_dim2 << 1) + 1] += ch1[(m2 + (
			ki + j * ch1_dim3) * ch1_dim2 << 1) + 1];
		cc1[(m1 + (ki + cc1_dim3) * cc1_dim2 << 1) + 2] += ch1[(m2 + (
			ki + j * ch1_dim3) * ch1_dim2 << 1) + 2];
/* L117: */
	    }
	}
/* L118: */
    }
    i__3 = ipph;
    for (l = 2; l <= i__3; ++l) {
	lc = ipp2 - l;
	i__2 = *lid;
	for (ki = 1; ki <= i__2; ++ki) {
	    m2 = m2s;
	    i__1 = m1d;
	    i__4 = *im1;
	    for (m1 = 1; i__4 < 0 ? m1 >= i__1 : m1 <= i__1; m1 += i__4) {
		m2 += *im2;
		cc1[(m1 + (ki + l * cc1_dim3) * cc1_dim2 << 1) + 1] = ch1[(m2 
			+ (ki + ch1_dim3) * ch1_dim2 << 1) + 1] + wa[(l - 1 + 
			wa_dim2) * wa_dim1 + 1] * ch1[(m2 + (ki + (ch1_dim3 <<
			 1)) * ch1_dim2 << 1) + 1];
		cc1[(m1 + (ki + lc * cc1_dim3) * cc1_dim2 << 1) + 1] = wa[(l 
			- 1 + (wa_dim2 << 1)) * wa_dim1 + 1] * ch1[(m2 + (ki 
			+ *ip * ch1_dim3) * ch1_dim2 << 1) + 1];
		cc1[(m1 + (ki + l * cc1_dim3) * cc1_dim2 << 1) + 2] = ch1[(m2 
			+ (ki + ch1_dim3) * ch1_dim2 << 1) + 2] + wa[(l - 1 + 
			wa_dim2) * wa_dim1 + 1] * ch1[(m2 + (ki + (ch1_dim3 <<
			 1)) * ch1_dim2 << 1) + 2];
		cc1[(m1 + (ki + lc * cc1_dim3) * cc1_dim2 << 1) + 2] = wa[(l 
			- 1 + (wa_dim2 << 1)) * wa_dim1 + 1] * ch1[(m2 + (ki 
			+ *ip * ch1_dim3) * ch1_dim2 << 1) + 2];
/* L113: */
	    }
	}
	i__4 = ipph;
	for (j = 3; j <= i__4; ++j) {
	    jc = ipp2 - j;
	    idlj = (l - 1) * (j - 1) % *ip;
	    war = wa[(idlj + wa_dim2) * wa_dim1 + 1];
	    wai = wa[(idlj + (wa_dim2 << 1)) * wa_dim1 + 1];
	    i__1 = *lid;
	    for (ki = 1; ki <= i__1; ++ki) {
		m2 = m2s;
		i__2 = m1d;
		i__5 = *im1;
		for (m1 = 1; i__5 < 0 ? m1 >= i__2 : m1 <= i__2; m1 += i__5) {
		    m2 += *im2;
		    cc1[(m1 + (ki + l * cc1_dim3) * cc1_dim2 << 1) + 1] += 
			    war * ch1[(m2 + (ki + j * ch1_dim3) * ch1_dim2 << 
			    1) + 1];
		    cc1[(m1 + (ki + lc * cc1_dim3) * cc1_dim2 << 1) + 1] += 
			    wai * ch1[(m2 + (ki + jc * ch1_dim3) * ch1_dim2 <<
			     1) + 1];
		    cc1[(m1 + (ki + l * cc1_dim3) * cc1_dim2 << 1) + 2] += 
			    war * ch1[(m2 + (ki + j * ch1_dim3) * ch1_dim2 << 
			    1) + 2];
		    cc1[(m1 + (ki + lc * cc1_dim3) * cc1_dim2 << 1) + 2] += 
			    wai * ch1[(m2 + (ki + jc * ch1_dim3) * ch1_dim2 <<
			     1) + 2];
/* L114: */
		}
	    }
/* L115: */
	}
/* L116: */
    }
    if (*ido > 1 || *na == 1) {
	goto L136;
    }
    i__3 = ipph;
    for (j = 2; j <= i__3; ++j) {
	jc = ipp2 - j;
	i__4 = *lid;
	for (ki = 1; ki <= i__4; ++ki) {
	    i__5 = m1d;
	    i__2 = *im1;
	    for (m1 = 1; i__2 < 0 ? m1 >= i__5 : m1 <= i__5; m1 += i__2) {
		chold1 = cc1[(m1 + (ki + j * cc1_dim3) * cc1_dim2 << 1) + 1] 
			- cc1[(m1 + (ki + jc * cc1_dim3) * cc1_dim2 << 1) + 2]
			;
		chold2 = cc1[(m1 + (ki + j * cc1_dim3) * cc1_dim2 << 1) + 1] 
			+ cc1[(m1 + (ki + jc * cc1_dim3) * cc1_dim2 << 1) + 2]
			;
		cc1[(m1 + (ki + j * cc1_dim3) * cc1_dim2 << 1) + 1] = chold1;
		cc1[(m1 + (ki + jc * cc1_dim3) * cc1_dim2 << 1) + 2] = cc1[(
			m1 + (ki + j * cc1_dim3) * cc1_dim2 << 1) + 2] - cc1[(
			m1 + (ki + jc * cc1_dim3) * cc1_dim2 << 1) + 1];
		cc1[(m1 + (ki + j * cc1_dim3) * cc1_dim2 << 1) + 2] += cc1[(
			m1 + (ki + jc * cc1_dim3) * cc1_dim2 << 1) + 1];
		cc1[(m1 + (ki + jc * cc1_dim3) * cc1_dim2 << 1) + 1] = chold2;
/* L119: */
	    }
	}
/* L120: */
    }
    return 0;
L136:
    i__3 = *lid;
    for (ki = 1; ki <= i__3; ++ki) {
	m2 = m2s;
	i__2 = m1d;
	i__5 = *im1;
	for (m1 = 1; i__5 < 0 ? m1 >= i__2 : m1 <= i__2; m1 += i__5) {
	    m2 += *im2;
	    ch1[(m2 + (ki + ch1_dim3) * ch1_dim2 << 1) + 1] = cc1[(m1 + (ki + 
		    cc1_dim3) * cc1_dim2 << 1) + 1];
	    ch1[(m2 + (ki + ch1_dim3) * ch1_dim2 << 1) + 2] = cc1[(m1 + (ki + 
		    cc1_dim3) * cc1_dim2 << 1) + 2];
/* L137: */
	}
    }
    i__5 = ipph;
    for (j = 2; j <= i__5; ++j) {
	jc = ipp2 - j;
	i__2 = *lid;
	for (ki = 1; ki <= i__2; ++ki) {
	    m2 = m2s;
	    i__3 = m1d;
	    i__4 = *im1;
	    for (m1 = 1; i__4 < 0 ? m1 >= i__3 : m1 <= i__3; m1 += i__4) {
		m2 += *im2;
		ch1[(m2 + (ki + j * ch1_dim3) * ch1_dim2 << 1) + 1] = cc1[(m1 
			+ (ki + j * cc1_dim3) * cc1_dim2 << 1) + 1] - cc1[(m1 
			+ (ki + jc * cc1_dim3) * cc1_dim2 << 1) + 2];
		ch1[(m2 + (ki + jc * ch1_dim3) * ch1_dim2 << 1) + 1] = cc1[(
			m1 + (ki + j * cc1_dim3) * cc1_dim2 << 1) + 1] + cc1[(
			m1 + (ki + jc * cc1_dim3) * cc1_dim2 << 1) + 2];
		ch1[(m2 + (ki + jc * ch1_dim3) * ch1_dim2 << 1) + 2] = cc1[(
			m1 + (ki + j * cc1_dim3) * cc1_dim2 << 1) + 2] - cc1[(
			m1 + (ki + jc * cc1_dim3) * cc1_dim2 << 1) + 1];
		ch1[(m2 + (ki + j * ch1_dim3) * ch1_dim2 << 1) + 2] = cc1[(m1 
			+ (ki + j * cc1_dim3) * cc1_dim2 << 1) + 2] + cc1[(m1 
			+ (ki + jc * cc1_dim3) * cc1_dim2 << 1) + 1];
/* L134: */
	    }
	}
/* L135: */
    }
    if (*ido == 1) {
	return 0;
    }
    i__5 = *ido;
    for (i__ = 1; i__ <= i__5; ++i__) {
	i__4 = *l1;
	for (k = 1; k <= i__4; ++k) {
	    m2 = m2s;
	    i__3 = m1d;
	    i__2 = *im1;
	    for (m1 = 1; i__2 < 0 ? m1 >= i__3 : m1 <= i__3; m1 += i__2) {
		m2 += *im2;
		cc[(m1 + (k + (i__ * cc_dim4 + 1) * cc_dim3) * cc_dim2 << 1) 
			+ 1] = ch[(m2 + (k + (i__ + ch_dim4) * ch_dim3) * 
			ch_dim2 << 1) + 1];
		cc[(m1 + (k + (i__ * cc_dim4 + 1) * cc_dim3) * cc_dim2 << 1) 
			+ 2] = ch[(m2 + (k + (i__ + ch_dim4) * ch_dim3) * 
			ch_dim2 << 1) + 2];
/* L130: */
	    }
	}
/* L131: */
    }
    i__5 = *ip;
    for (j = 2; j <= i__5; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    m2 = m2s;
	    i__3 = m1d;
	    i__4 = *im1;
	    for (m1 = 1; i__4 < 0 ? m1 >= i__3 : m1 <= i__3; m1 += i__4) {
		m2 += *im2;
		cc[(m1 + (k + (j + cc_dim4) * cc_dim3) * cc_dim2 << 1) + 1] = 
			ch[(m2 + (k + (j * ch_dim4 + 1) * ch_dim3) * ch_dim2 
			<< 1) + 1];
		cc[(m1 + (k + (j + cc_dim4) * cc_dim3) * cc_dim2 << 1) + 2] = 
			ch[(m2 + (k + (j * ch_dim4 + 1) * ch_dim3) * ch_dim2 
			<< 1) + 2];
/* L122: */
	    }
	}
/* L123: */
    }
    i__5 = *ip;
    for (j = 2; j <= i__5; ++j) {
	i__4 = *ido;
	for (i__ = 2; i__ <= i__4; ++i__) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		m2 = m2s;
		i__2 = m1d;
		i__1 = *im1;
		for (m1 = 1; i__1 < 0 ? m1 >= i__2 : m1 <= i__2; m1 += i__1) {
		    m2 += *im2;
		    cc[(m1 + (k + (j + i__ * cc_dim4) * cc_dim3) * cc_dim2 << 
			    1) + 1] = wa[i__ + (j - 1 + wa_dim2) * wa_dim1] * 
			    ch[(m2 + (k + (i__ + j * ch_dim4) * ch_dim3) * 
			    ch_dim2 << 1) + 1] - wa[i__ + (j - 1 + (wa_dim2 <<
			     1)) * wa_dim1] * ch[(m2 + (k + (i__ + j * 
			    ch_dim4) * ch_dim3) * ch_dim2 << 1) + 2];
		    cc[(m1 + (k + (j + i__ * cc_dim4) * cc_dim3) * cc_dim2 << 
			    1) + 2] = wa[i__ + (j - 1 + wa_dim2) * wa_dim1] * 
			    ch[(m2 + (k + (i__ + j * ch_dim4) * ch_dim3) * 
			    ch_dim2 << 1) + 2] + wa[i__ + (j - 1 + (wa_dim2 <<
			     1)) * wa_dim1] * ch[(m2 + (k + (i__ + j * 
			    ch_dim4) * ch_dim3) * ch_dim2 << 1) + 1];
/* L124: */
		}
	    }
/* L125: */
	}
/* L126: */
    }
    return 0;
} /* cmfgkb_ */

