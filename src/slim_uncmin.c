/* The nonlinear minimizer is written after the source code of The Dennis + Schnabel Minimizer -- used by R's  nlm() 
 * 
 R : A Computer Language for Statistical Data Analysis
 *  R is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h> /* DBL_MAX */
#include <cblas.h> /* ddot, dnrm2, dscal */
#include "arrays.h"
#include "uncmin.h"

///////////////////////////////////////////////////////////////////////
#ifndef FALSE
#define  FALSE 0
#endif
#ifndef TRUE
#define  TRUE 1
#endif

typedef int bool;


#define Rexp10(x) pow(10.0, x)
////////////////////////////////////////////////////////////////////////
static double FMAX2(double x, double y)
{
	return (x < y) ? y : x;
}

static double FMIN2(double x, double y)
{
	return (x < y) ? x : y;
}


///////////////////////////////////////////////

static void lltslv(int nr, int n, double *a, double *x, double *b)
{
/*	solve ax=b where a has the form l(l-transpose)
          LL^Tx=b
 *	but only the lower triangular part, l, is stored.
 * PARAMETERS :
 *	nr	     --> row dimension of matrix
 *	n	     --> dimension of problem
 *	a(n,n)	     --> matrix of form l(l-transpose).
 *			 on return a is unchanged.
 *	x(n)	    <--	 solution vector
 *	b(n)	     --> right-hand side vector
 *	note
 *	if b is not required by calling program, then
 *	b and x may share the same storage. */
	if (x != b) cblas_dcopy (n, b, 1, x,1);
    cblas_dtrsv (CblasColMajor, CblasLower , CblasNoTrans, CblasNonUnit, n, a, nr,x,1);
    cblas_dtrsv (CblasColMajor, CblasLower , CblasTrans,   CblasNonUnit, n, a, nr,x,1);

}

///////////////////////////////////////////////

static void qraux1(int nr, int n, double *r, int i)
{
/* Interchange rows i,i+1 of the upper hessenberg matrix r, columns i to n .
 * PARAMETERS :
 *	nr	     --> row dimension of matrix
 *	n	     --> dimension of matrix
 *	r[n*n]	    <--> upper hessenberg matrix
 *	i	     --> index of row to interchange (i < n-1)
*/
  double tmp;
  double *r1, *r2;
  /* pointer arithmetic : */
  r1 = r + i + i * nr;
  r2 = r1 + 1;

  while(n-- > i) {
    tmp = *r1; *r1 = *r2; *r2 = tmp;
    r1 += nr;
    r2 += nr;
  }
} /* qraux1 */

static void qraux2(int nr, int n, double *r, int i, double a, double b)
{
/* Pre-multiply r by the jacobi rotation j(i,i+1,a,b) .
 * PARAMETERS :
 *	nr	     --> row dimension of matrix
 *	n	     --> dimension of matrix
 *	r(n,n)	    <--> upper hessenberg matrix
 *	i	     --> index of row
 *	a	     --> scalar
 *	b	     --> scalar */

  double c, s;
  double y, z, den;
  double *r1, *r2;

  den = hypot(a,b);
  c = a / den;
  s = b / den;

  /* pointer arithmetic : */
  r1 = r + i + i*nr;
  r2 = r1 + 1;

  while(n-- > i) {
    y = *r1;
    z = *r2;
    *r1 = c * y - s * z;
    *r2 = s * y + c * z;
    r1 += nr;
    r2 += nr;
  }
} /* qraux2 */

static void
qrupdt(int nr, int n, double *a, double *u, double *v)
{
/* Find an orthogonal (n*n) matrix (q*) and an upper triangular (n*n)
 * matrix (r*) such that (q*)(r*)=r+u(v+)
 * PARAMETERS :
 *	nr	     --> row dimension of matrix
 *	n	     --> dimension of problem
 *	a(n,n)	    <--> on input:  contains r
 *			 on output: contains (r*)
 *	u(n)	     --> vector
 *	v(n)	     --> vector */

    int i, j, k;
    double t1, t2;
    int ii;

    /*	determine last non-zero in u(.) */

    for(k = n-1; k > 0 && u[k] == 0.0; k--)
	;

    /*	(k-1) jacobi rotations transform
     *	    r + u(v+) --> (r*) + (u(1)*e1)(v+)
     *	which is upper hessenberg */

    if (k > 0)
    {
    	ii = k;
    	while(ii > 0)
    	{
    		i = ii - 1;
    		if (u[i] == 0.0)
    		{
    			qraux1(nr, n, a, i);
    			u[i] = u[ii];
    		}
    		else
    		{
    			qraux2(nr, n, a, i, u[i], -u[ii]);
    			u[i] = hypot(u[i], u[ii]);
    		}
    		ii = i;
    	}
    }
    /*	r <-- r + (u(1)*e1)(v+) */
    for (j = 0; j < n; ++j)	a[j * nr] += u[0] * v[j];
    /*	(k-1) jacobi rotations transform upper hessenberg r
     *	to upper triangular (r*) */

    for (i = 0; i < k; ++i)
    {
    	if (a[i + i * nr] == 0.)   qraux1(nr, n, a, i);
    	else
    	{
    		t1 = a[i + i * nr];
    		t2 = -a[i + 1 + i * nr];
    		qraux2(nr, n, a, i, t1, t2);
    	}
    }
} /* qrupdt */


static void
lnsrch(int n, double *x, double f, double *g, double *p, double *xpls,
       double *fpls, fcn_p fcn, void *state, bool *mxtake, int *iretcd,
       double stepmx, double steptl)
{
/* A.6.3.1 p. 325
 * Find a next newton iterate by line search.  (iff  method == 1)
 * PARAMETERS :
 *	n	     --> dimension of problem
 *	x(n)	     --> old iterate:	x[k-1]
 *	f	     --> function value at old iterate, f(x)
 *	g(n)	     --> gradient at old iterate, g(x), or approximate
 *	p(n)	     --> non-zero newton step
 *	xpls(n)	    <--	 new iterate x[k]
 *	fpls	    <--	 function value at new iterate, f(xpls)
 *	fcn	     --> name of subroutine to evaluate function
 *	state	    <--> information other than x and n that fcn requires.
 *			 state is not modified in lnsrch (but can be
 *			 modified by fcn).
 *	iretcd	    <--	 return code
 *	mxtake	    <--  boolean flag indicating step of maximum length used
 *	stepmx	     --> maximum allowable step size
 *	steptl	     --> relative step size at which successive iterates
 *			 considered close enough to terminate algorithm
 *	internal variables
 *	sln		 newton length
 *	rln		 relative length of newton step
*/

    int i, one = 1;
    bool firstback = TRUE;
    double disc;
    double a3, b;
    double t1, t2, t3, lambda, tlmbda, rmnlmb;
    double scl, rln, sln, slp;
    double temp1;
    double pfpls = 0., plmbda = 0.; /* -Wall */

    temp1 = 0.;
    for (i = 0; i < n; ++i)
	temp1 +=  p[i] * p[i];
    sln = sqrt(temp1);
    if (sln > stepmx) {
	/*	newton step longer than maximum allowed */
	scl = stepmx / sln;
	cblas_dscal(n, scl, p, one);
	sln = stepmx;
    }
	slp = cblas_ddot(n, g, one, p, one);
    rln = 0.;
    for (i = 0; i < n; ++i) {
	temp1 = fabs(p[i])/ FMAX2(fabs(x[i]), 1);
	if(rln < temp1) rln = temp1;
    }
    rmnlmb = steptl / rln;
    lambda = 1.0;

    /*	check if new iterate satisfactory.  generate new lambda if necessary. */

    *mxtake = FALSE;
    *iretcd = 2;
    do 
	{
		for (i = 0; i < n; ++i)	    xpls[i] = x[i] + lambda * p[i];
		(*fcn)(n, xpls, fpls, state);
		if (*fpls <= f + slp * 1e-4 * lambda) 
		{ /* solution found */
			*iretcd = 0;
			if (lambda == 1. && sln > stepmx * .99) *mxtake = TRUE;
			return;
		}
		/* else : solution not (yet) found */
		/* First find a point with a finite value */

		if (lambda < rmnlmb) 
		{
			/* no satisfactory xpls found sufficiently distinct from x */
			*iretcd = 1;
			return;
		}
		else 
		{ 	/*	calculate new lambda */
			/* modifications by BDR 2000/01/05 to cover non-finite values
			* ">=" instead of "==" :  MM 2001/07/24 */
			if (*fpls >= DBL_MAX) 
			{
				lambda *= 0.1;
				firstback = TRUE;
			}
			else 
			{
				if (firstback) 
				{ /*	first backtrack: quadratic fit */
					tlmbda = -lambda * slp / ((*fpls - f - slp) * 2.);
					firstback = FALSE;
				}
				else 
				{ /*	all subsequent backtracks: cubic fit */
					t1 = *fpls - f - lambda * slp;
					t2 = pfpls - f - plmbda * slp;
					t3 = 1. / (lambda - plmbda);
					a3 = 3. * t3 * (t1 / (lambda * lambda) - t2 / (plmbda * plmbda));
					b = t3 * (t2 * lambda / (plmbda * plmbda) - t1 * plmbda / (lambda * lambda));
					disc = b * b - a3 * slp;
					if (disc > b * b)
					/* only one positive critical point, must be minimum */
						tlmbda = (-b + ((a3 < 0)? -sqrt(disc): sqrt(disc))) /a3;
					else
						/* both critical points positive, first is minimum */
						tlmbda = (-b + ((a3 < 0)? sqrt(disc): -sqrt(disc))) /a3;

					if (tlmbda > lambda * .5)		tlmbda = lambda * .5;
				}
				plmbda = lambda;
				pfpls = *fpls;
				if (tlmbda < lambda * .1)    lambda *= .1;
				else	    lambda = tlmbda;
			}
		}
    } while(*iretcd > 1);
} /* lnsrch */



static void
secfac(int nr, int n, double *x, double *g, double *a, double *xpls,
       double *gpls, double epsm, int itncnt, double rnf, int iagflg,
       bool *noupdt, double *s, double *y, double *u, double *w)
{
/* A9.4.2
 * Update hessian by the bfgs factored method  (only when  method == 1 or 2)
 * PARAMETERS :
 *	nr	     --> row dimension of matrix
 *	n	     --> dimension of problem
 *	x(n)	     --> old iterate, x[k-1]
 *	g(n)	     --> gradient or approximate at old iterate
 *	a(n,n)	    <--> on entry: cholesky decomposition of hessian in
 *			   lower part and diagonal.
 *			 on exit:  updated cholesky decomposition of hessian
 *			   in lower triangular part and diagonal
 *	xpls(n)	     --> new iterate, x[k]
 *	gpls(n)	     --> gradient or approximate at new iterate
 *	epsm	     --> machine epsilon
 *	itncnt	     --> iteration count
 *	rnf	     --> relative noise in optimization function fcn
 *	iagflg	     --> =1 if analytic gradient supplied, =0 itherwise
 *	noupdt	    <--> boolean: no update yet
 *			 [retain value between successive calls]
 *	s(n)	     --> workspace
 *	y(n)	     --> workspace
 *	u(n)	     --> workspace
 *	w(n)	     --> workspace */

    double ynrm2;
    int i, j, jnr;
    int skpupd;
    double snorm2, reltol;
    double alp, den1, den2;
    *noupdt = (itncnt == 1);
    for (i = 0; i < n; ++i)
    {
    	s[i] = xpls[i] - x[i];
    	y[i] = gpls[i] - g[i];
    }
	den1   = cblas_ddot(n, s, 1, y, 1);
    snorm2 = cblas_dnrm2(n, s, 1);
    ynrm2  = cblas_dnrm2(n, y, 1);

    if (den1 < sqrt(epsm) * snorm2 * ynrm2)
	return;

   //mvmltu(nr, n, a, s, u);
   cblas_dcopy (n, s, 1, u,1);
   cblas_dtrmv (CblasColMajor, CblasLower , CblasTrans, CblasNonUnit, n, a, nr, u, 1);
	den2 = cblas_ddot(n, u, 1, u, 1);

    /*	l <-- sqrt(den1/den2)*l */

    alp = sqrt(den1 / den2);
    if (*noupdt)
    {
    	for (j = 0; j < n; ++j)
    	{
    		u[j] = alp * u[j];
	    	for (i = j; i < n; ++i)
	    	{
	    		jnr= j * nr;
	    		a[i + jnr] *= alp;
	    	}
    	}
    	*noupdt = FALSE;
    	den2 = den1;
    	alp = 1.;
    }
    /*	w = l(l+)s = hs */
    //mvmltl(nr, n, a, u, w);
    cblas_dcopy (n, u, 1, w,1);
    cblas_dtrmv (CblasColMajor, CblasLower , CblasNoTrans, CblasNonUnit, n, a, nr, w, 1);
    if (iagflg == 0)
	reltol = sqrt(rnf);
    else
	reltol = rnf;

    skpupd = TRUE;
    for (i = 0; i < n; ++i)
    {
    	skpupd = (fabs(y[i] - w[i]) <  reltol * FMAX2(fabs(g[i]), fabs(gpls[i])));
    	if(!skpupd)	    break;
    }
    if(skpupd)	return;

    /*	  w = y-alp*l(l+)s */
    for (i = 0; i < n; ++i)
	w[i] = y[i] - alp * w[i];


    /*	  alp=1/sqrt(den1*den2) */
    alp /= den1;

    /*	  u=(l+)/sqrt(den1*den2) = (l+)s/sqrt((y+)s * (s+)l(l+)s) */

    for (i = 0; i < n; ++i)
	u[i] *= alp;
    /*	  copy l into upper triangular part.  zero l. */
    int inr, ijnr;
    for (i = 1; i < n; ++i)
    {
    	for (j = 0; j < i; ++j)
    	{
    		inr= i * nr;
    		ijnr= i+ j * nr;
    		a[j + inr] = a[ijnr];
    		a[ijnr] = 0.;
    	}
    }

    /*	  find q, (l+) such that  q(l+) = (l+) + u(w+) */

    qrupdt(nr, n, a, u, w);

    /* upper triangular part and diagonal of a[] now contain updated
     * cholesky decomposition of hessian.
     * copy back to lower triangular part.
     */
    for (i = 1; i < n; ++i)
	for (j = 0; j < i; ++j)	    a[i + j * nr] = a[j + i * nr];
} /* secfac */


static void
fstofd(int n, double *xpls, fcn_p fcn, void *state,  double fpls, double *g, double rnoise)
{
/*	find first derivative (gradient) forward finite difference approximation "g" to the
 *	first derivative of the function fcn  evaluated at the new iterate "xpls".
 * PARAMETERS :
 *	n	     --> number of columns in a; dimension of problem
 *	xpls(n)	     --> new iterate:  x[k]
 *	fcn	     --> name of subroutine to evaluate function
 *	state	    <--> information other than x and n that fcn requires.
 *			 state is not modified in fstofd (but can be
 *			 modified by fcn).
 *	g(n)	    <--	 finite difference approximation (see note).  only
 *			 lower triangular matrix and diagonal are returned
  *	rnoise	     --> relative noise in fcn [f(x)]
 *	internal variables
 *	stepsz - stepsize in the j-th variable direction
 */

    int i, j,jnr;
    double fhat, xtmpj, stepsz;
    double sroot_rnoise= sqrt(rnoise);
    for (j = 0; j < n; ++j)
    {
    	stepsz = sroot_rnoise * FMAX2(fabs(xpls[j]), 1);
    	xtmpj = xpls[j];
    	xpls[j] = xtmpj + stepsz;
    	(*fcn)(n, xpls, &fhat, state);
    	xpls[j] = xtmpj;
   		g[j] = (fhat - fpls) / stepsz;
    }
 } /* fstofd */
//////////////////////////////////////////////////////////////////////////////////////
static void fstocd(int n, double *x, fcn_p fcn, void *state,
		    double rnoise, double *g)
{
/* Find central difference approximation g to the first derivative
 * (gradient) of the function defined by fcn at the point x.
 * PARAMETERS :
 *	n	     --> dimension of problem
 *	x	     --> point at which gradient is to be approximated.
 *	fcn	     --> name of subroutine to evaluate function.
 *	state	    <--> information other than x and n that fcn requires.
 *			 state is not modified in fstocd (but can be
 *			 modified by fcn).
 *	rnoise	     --> relative noise in fcn [f(x)].
 *	g	    <--	 central difference approximation to gradient.
 */
    int i;
    double stepi, fplus, fminus, xtempi;

    /*	find i th  stepsize, evaluate two neighbors in direction of i th */
    /*	unit vector, and evaluate i th	component of gradient. */
    double ttt= pow(rnoise, 1.0/3.0);
    for (i = 0; i < n; ++i) {
	xtempi = x[i];
	stepi = ttt * FMAX2(fabs(xtempi), 1);
	x[i] = xtempi + stepi;
	(*fcn)(n, x, &fplus, state);
	x[i] = xtempi - stepi;
	(*fcn)(n, x, &fminus, state);
	x[i] = xtempi;
	g[i] = (fplus - fminus) / (stepi * 2.);
    }
} /* fstocd */


static
int opt_stop(int n, double *xpls, double fpls, double *gpls, double *x,
	     int itncnt, int *icscmx, double gradtl, double steptl,
	       int itnlim,
	     int iretcd, int mxtake, int *msg)
{
/* Unconstrained minimization stopping criteria :
 * Find whether the algorithm should terminate, due to any
 * of the following:
 *	1) problem solved within user tolerance
 *	2) convergence within user tolerance
 *	3) iteration limit reached
 *	4) divergence or too restrictive maximum step (stepmx) suspected
 * ARGUMENTS :
 *	n	     --> dimension of problem
 *	xpls(n)	     --> new iterate x[k]
 *	fpls	     --> function value at new iterate f(xpls)
 *	gpls(n)	     --> gradient at new iterate, g(xpls), or approximate
 *	x(n)	     --> old iterate x[k-1]
 *	itncnt	     --> current iteration k
 *	icscmx	    <--> number consecutive steps >= stepmx
 *			 [retain value between successive calls]
 *	gradtl	     --> tolerance at which relative gradient considered close
 *			 enough to zero to terminate algorithm
 *	steptl	     --> relative step size at which successive iterates
 *			 considered close enough to terminate algorithm
 *	sx(n)	     --> diagonal scaling matrix for x
 *	fscale	     --> estimate of scale of objective function
 *	itnlim	     --> maximum number of allowable iterations
 *	iretcd	     --> return code
 *	mxtake	     --> boolean flag indicating step of maximum length used
 *	msg	     --> if msg includes a term 8, suppress output
 *
 * VALUE :
 *	`itrmcd' : termination code
 */

    int i, jtrmcd;
    double d, relgrd, relstp, rgx, rsx;

    /*	last global step failed to locate a point lower than x */
    if (iretcd == 1)
	return 3;

    /* else : */

    /* find direction in which relative gradient maximum. */

    /* check whether within tolerance */
    d = FMAX2(fabs(fpls), 1);
    rgx = 0.;
    for (i = 0; i < n; ++i) 
	{
		relgrd = fabs(gpls[i]) * FMAX2(fabs(xpls[i]), 1) / d;
		if(rgx < relgrd) rgx = relgrd;
    }
    jtrmcd = 1;

    if (rgx > gradtl) 
	{


		if (itncnt == 0)
			return 0;

		/* find direction in which relative stepsize maximum */

		/* check whether within tolerance. */
		rsx = 0.;

		for (i = 0; i < n; ++i) 
		{
			relstp = fabs(xpls[i] - x[i]) / FMAX2(fabs(xpls[i]), 1.0);
			if(rsx < relstp) rsx = relstp;
		}
//printf("relstp .... steptl: %g  %g \n", rsx,  steptl );
		jtrmcd = 2;
		if (rsx > steptl) 
		{ /*	check iteration limit */
			jtrmcd = 4;
			if (itncnt < itnlim)
			{
				/*	check number of consecutive steps \ stepmx */
				if (!mxtake) 
				{
					*icscmx = 0; return 0;
				} 
				else 
				{
					++(*icscmx);
					if (*icscmx < 5) return 0;
					jtrmcd = 5;
				}
			}
		}
    }

    return jtrmcd;
} /* opt_stop */


static void
prt_result( int n, const double x[], double f, const double g[],
	   const double *a, const double p[], int itncnt, int iflg)
{
    printf("iteration = %d\n", itncnt);
    printf("Parameter:\n");
	int i;
	for (i = 0; i < n; ++i) printf("%f ", x[i]);
	printf("\n");
     printf("Function Value: %f \n", f );
  	printf("\n");
} /* prt_result */

static void
optdrv_end( int n, double *xpls, double *x, double *gpls,
	   double *g, double *fpls, double f, double *a, double *p,
	   int itncnt, int itrmcd, int *msg)
{
    int i;

    /*	termination :
	reset xpls,fpls,gpls,  if previous iterate solution */
    if (itrmcd == 3)
    {
    	*fpls = f;
    	for (i = 0; i < n; ++i)
    	{
    		xpls[i] = x[i];
    		gpls[i] = g[i];
    	}
    }
 //   if (*msg / 8 % 2 == 0)
//	(*print_result)(nr, n, xpls, *fpls, gpls, a, p, itncnt, 0);

    *msg = 0;
} /* optdrv_end */


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
static void printmes(char *message)
{
	printf("Error: %s\n",message);
}

static void opterror(int nerr)
/* Fatal errors - we don't deliver an answer */
{
    switch(nerr) 
	{
    case -1:
		printmes("non-positive number of parameters in nlm");
		break;
    case -2:
		printmes("nlm is inefficient for 1-d problems");
		break;
    case -3:
		printmes("invalid gradient tolerance in nlm");
		break;
    case -4:
		printmes("invalid iteration limit in nlm");
		break;
    case -5:
		printmes("minimization function has no good digits in nlm");
		break;
    default:
		printmes("*** unknown error message  in nlm()\n*** should not happen!");
    }
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
	

static void optcode(int code)
/* Warnings - we return a value, but print a warning */
{
    switch(code) {
    case 1:
		printf("Relative gradient close to zero.\n");
		printf("Current iterate is probably solution.\n");
		break;
    case 2:
		printf("Successive iterates within tolerance.\n");
		printf("Current iterate is probably solution.\n");
		break;
    case 3:
		printf("Last global step failed to locate a point lower than x.\n");
		printf("Either x is an approximate local minimum of the function,\n\
				the function is too non-linear for this algorithm,\n\
				or steptol is too large.\n");
		break;
    case 4:
		printf("Iteration limit exceeded.  Algorithm failed.\n");
		break;
    case 5:
		printf("Maximum step size exceeded 5 consecutive times.\n\
				Either the function is unbounded below,\n\
				becomes asymptotic to a finite value\n\
				from above in some direction,\n"\
				"or stepmx is too small.\n");
		break;
    }
    printf("\n");
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void nlmin_slim( fcn_p fcn , int n , double *x, void *state, double *fpls ,int itnlim ,  int verbose )
/****************************************************************************************************
fcn_p fcn   <-->  type of pointer to the target and gradient functions it is should be delcared as 
                  void fcn(int n, double *x, double *fpls, void *)
x(n)	    <--> on exit:  x should contain the local minimum
fpls	    <--> on exit:  function value at solution (the local minimum)	 
itnlim	     --> maximum number of allowable iterations
verbose = 1 print the warnnnig meessege
*****************************************************************************************************/

{
	double   f,rnf, gradtl, stepmx, steptl,  dlt;
    int itrmcd, i, j, k,  msg, ndigit,  itncnt;
	int iagflg=0;
	bool mxtake = FALSE, noupdt;
    int  iretcd = 0, icscmx = 0;
    itncnt = 0;
	//////// set tolerances ///////
	double epsm = DBL_EPSILON;		////// for IEEE : = 2^-52	  ~= 2.22  e-16
	gradtl = pow(epsm, 1./3.);		////// for IEEE : = 2^(-52/3) ~= 6.055 e-6
	steptl = sqrt(epsm);		    ////// for IEEE : = 2^-26	  ~= 1.490 e-8
	//steptol=1e-6;
	////// compute default maximum step size if not provided
	double stpsiz = 0.;
	//for (i = 0; i < n; ++i)	    stpsiz += x[i] * x[i] ;
	//stepmx = 1000. * FMAX2(sqrt(stpsiz), 1.);
	stpsiz=cblas_dnrm2 (n, x, 1);
	stepmx = 1000. * FMAX2(stpsiz, 1.);
	ndigit = (int) (-log10(epsm));
    rnf = Rexp10(-(double)ndigit);
    rnf = FMAX2(rnf, epsm);
	//msg=9;
	//if (((msg/4) % 2) && !iahflg)  msg -= 4; /* skip check of analytic Hessian */
    //if (((msg/2) % 2) && !iagflg)   msg -= 2;  /* skip check of analytic gradient */
	msg=3;
  	ONE_DIM_ARRAY(xpls, double,n );
	ONE_DIM_ARRAY(gpls, double,n );
	ONE_DIM_ARRAY(a, double,(n*n) );
	ONE_DIM_ARRAY(wrk, double,(6*n) );
	double *g		= wrk ;
	double *p		= wrk + n ;
	double *wrk0	= wrk + n * 2;
	double *wrk1	= wrk + n * 3;
	double *wrk2	= wrk + n * 4;
	double *wrk3	= wrk + n * 5;
    for (i = 0; i < n; ++i)	p[i] = 0.;
    /////	evaluate fcn(x)
    (*fcn)(n, x, &f, state);
    /////evaluate finite difference gradient
	fstofd( n, x, (fcn_p)fcn, state, f, g, rnf);
	iretcd = -1;
    itrmcd = opt_stop(n, x, f, g, wrk1, itncnt, &icscmx, gradtl, steptl, itnlim, iretcd,FALSE, &msg);
    if (itrmcd != 0)
    {
    	optdrv_end(n, xpls, x, gpls, g, fpls, f, a, p, itncnt, 3, &msg);
    	return;
    }
	///// initialize hessian to be obtained by secant because optimization function expensive to evaluate
    for (i = 0; i < n; ++i)
     {
    	   a[i + i * n] = 1;
    	    for (j = 0; j < i; ++j)    a[i + j * n] = 0.;
     }
    //if (msg / 8 % 2 == 0)	prt_result( n, x, f, g, a, p, *itncnt, 1);
    //////////////////////// THE Iterations /////////////////////
    while(1)
    {
    	++(itncnt);
    	iagflg = 0;
    	L105:
		/////// solve for newton step: ap=-g ///////////
		for (i = 0; i < n; ++i) wrk1[i] = - g[i];
		lltslv(n, n, a, p, wrk1);
		lnsrch(n, x, f, g, p, xpls, fpls, (fcn_p)fcn, state, &mxtake, &iretcd, stepmx, steptl);
		//////	retry central difference gradient if could not find satisfactory step with forward difference gradient//////
		if (iretcd == 1 && iagflg == 0)
		{
			//printf("Central %d\n",itncnt);
			iagflg = -1; /////	 set iagflg for central differences ////
			fstocd(n, x, (fcn_p)fcn, state, rnf, g);
			goto L105;
		}
		/////	calculate step for output ///
		for (i = 0; i < n; ++i)	    p[i] = xpls[i] - x[i];
		////	calculate gradient at xpls /////
		if(iagflg==-1)
		{
			///// central difference gradient /////
			fstocd(n, xpls, (fcn_p)fcn, state,  rnf, gpls);
		}
		else
		{
			///// forward difference gradient/////
			fstofd( n, xpls, (fcn_p)fcn, state, *fpls, gpls,  rnf);
		}
		//////	check whether stopping criteria satisfied /////
		itrmcd = opt_stop(n, xpls, *fpls, gpls, x, itncnt, &icscmx,gradtl, steptl,   itnlim, iretcd, mxtake, &msg);
		if(itrmcd != 0) break;
		/////	evaluate hessian at xpls  expensive obj.fun.
		secfac(n, n, x, g, a, xpls, gpls, epsm, itncnt, rnf, iagflg, &noupdt, wrk0, wrk1, wrk2, wrk3);
		if (verbose >1 && msg / 16 % 2 == 1)  prt_result( n, xpls, *fpls, gpls, a, p, itncnt, 1);
		f = *fpls;

		for (i = 0; i < n; ++i)
		{
			x[i] = xpls[i];
			g[i] = gpls[i];
		}
    } /* END while(1) */
    optdrv_end( n, xpls, x, gpls, g, fpls, f, a, p, itncnt, itrmcd, &msg);


////////////////////////////////////////////////
	free(wrk);
	free(a);
	free(gpls);
    if (msg < 0)	opterror(msg);
    if (itrmcd != 0  && verbose>1)	optcode(itrmcd);
	for (i = 0; i < n; ++i)  x[i]=xpls[i];
	free(xpls);

	   
}
///////////////////////////////////////////////////////////////






