package org.rfoundation.R.nmath;

import org.rfoundation.R.Rmath;

public class Norm {
	/*
	 *  Mathlib : A C Library of Special Functions
	 *  Copyright (C) 1998 Ross Ihaka
	 *  Copyright (C) 2000	    The R Development Core Team
	 *  Copyright (C) 2003	    The R Foundation
	 *
	 *  This program is free software; you can redistribute it and/or modify
	 *  it under the terms of the GNU General Public License as published by
	 *  the Free Software Foundation; either version 2 of the License, or
	 *  (at your option) any later version.
	 *
	 *  This program is distributed in the hope that it will be useful,
	 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
	 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	 *  GNU General Public License for more details.
	 *
	 *  You should have received a copy of the GNU General Public License
	 *  along with this program; if not, a copy is available at
	 *  http://www.r-project.org/Licenses/
	 *
	 *  SYNOPSIS
	 *
	 *	double dnorm4(double x, double mu, double sigma, int give_log)
	 *	      {dnorm (..) is synonymous and preferred inside R}
	 *
	 *  DESCRIPTION
	 *
	 *	Compute the density of the normal distribution.
	 */
	public static double dnorm(double x, double mu, double sigma, boolean give_log) {
		return dnorm4(x, mu, sigma, give_log);
	}

	private static double dnorm4(double x, double mu, double sigma, boolean give_log) {
	    if(Double.isInfinite(sigma)) return DPQ.RD0(give_log);
	    if(Double.isInfinite(x) && mu == x) return Double.NaN;/* x-mu is NaN */
	    if (sigma <= 0) {
		if (sigma < 0) return Double.NaN;
		/* sigma == 0 */
		return (x == mu) ? Double.POSITIVE_INFINITY : DPQ.RD0(give_log);
	    }
	    x = (x - mu) / sigma;

	    if(Double.isInfinite(x)) return DPQ.RD0(give_log);
	    return (give_log ?
		    -(Rmath.M_LN_SQRT_2PI  +	0.5 * x * x + Math.log(sigma)) :
		    Rmath.M_1_SQRT_2PI * Math.exp(-0.5 * x * x)  /	  sigma);
	    /* M_1_SQRT_2PI = 1 / sqrt(2 * pi) */
	}

	/*
	 *  Mathlib : A C Library of Special Functions
	 *  Copyright (C) 1998	    Ross Ihaka
	 *  Copyright (C) 2000-2010 The R Development Core Team
	 *  Copyright (C) 2003	    The R Foundation
	 *
	 *  This program is free software; you can redistribute it and/or modify
	 *  it under the terms of the GNU General Public License as published by
	 *  the Free Software Foundation; either version 2 of the License, or
	 *  (at your option) any later version.
	 *
	 *  This program is distributed in the hope that it will be useful,
	 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
	 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	 *  GNU General Public License for more details.
	 *
	 *  You should have received a copy of the GNU General Public License
	 *  along with this program; if not, a copy is available at
	 *  http://www.r-project.org/Licenses/
	 *
	 *  SYNOPSIS
	 *
	 *   #include <Rmath.h>
	 *
	 *   double pnorm5(double x, double mu, double sigma, int lower_tail,int log_p);
	 *	   {pnorm (..) is synonymous and preferred inside R}
	 *
	 *   void   pnorm_both(double x, double *cum, double *ccum,
	 *		       int i_tail, int log_p);
	 *
	 *  DESCRIPTION
	 *
	 *	The main computation evaluates near-minimax approximations derived
	 *	from those in "Rational Chebyshev approximations for the error
	 *	function" by W. J. Cody, Math. Comp., 1969, 631-637.  This
	 *	transportable program uses rational functions that theoretically
	 *	approximate the normal distribution function to at least 18
	 *	significant decimal digits.  The accuracy achieved depends on the
	 *	arithmetic system, the compiler, the intrinsic functions, and
	 *	proper selection of the machine-dependent constants.
	 *
	 *  REFERENCE
	 *
	 *	Cody, W. D. (1993).
	 *	ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of
	 *	Special Function Routines and Test Drivers".
	 *	ACM Transactions on Mathematical Software. 19, 22-32.
	 *
	 *  EXTENSIONS
	 *
	 *  The "_both" , lower, upper, and log_p  variants were added by
	 *  Martin Maechler, Jan.2000;
	 *  as well as log1p() and similar improvements later on.
	 *
	 *  James M. Rath contributed bug report PR#699 and patches correcting SIXTEN
	 *  and if() clauses {with a bug: "|| instead of &&" -> PR #2883) more in line
	 *  with the original Cody code.
	 */
	/**
	 * Wrapper class that enables passing by reference in Java
	 * @author Yogesh Saletore
	 *
	 * @param <E>
	 */
	private static class Wrapper<E> {
		E obj;
		
		public Wrapper(E o) {
			obj = o;
		}
	}
	
	public static double pnorm(double x, double mu, double sigma, boolean lower_tail, boolean log_p) {
		return pnorm5(x, mu, sigma, lower_tail, log_p);
	}

	private static double pnorm5(double x, double mu, double sigma, boolean lower_tail, boolean log_p) {
	    double p = 0, cp = 0;

	    /* Note: The structure of these checks has been carefully thought through.
	     * For example, if x == mu and sigma == 0, we get the correct answer 1.
	     */
	    if(Double.isInfinite(x) && mu == x) return Double.NaN;/* x-mu is NaN */
	    if (sigma <= 0) {
		if(sigma < 0) return Double.NaN;
		/* sigma = 0 : */
		return (x < mu) ? DPQ.R_DT_0(lower_tail, log_p) : DPQ.R_DT_1(lower_tail, log_p);
	    }
	    p = (x - mu) / sigma;
	    if(Double.isInfinite(p))
		return (x < mu) ? DPQ.R_DT_0(lower_tail, log_p) : DPQ.R_DT_1(lower_tail, log_p);
	    x = p;

	    Wrapper<Double> pWrap = new Wrapper<Double>(p);
	    Wrapper<Double> cpWrap = new Wrapper<Double>(cp);
	    pnorm_both(x, pWrap, cpWrap, !lower_tail, log_p);

	    return(lower_tail ? pWrap.obj : cpWrap.obj);
	}

	private static final int SIXTEN = 16; /* Cutoff allowing exact "*" and "/" */
	
	private static final double[] a = {
		2.2352520354606839287,
		161.02823106855587881,
		1067.6894854603709582,
		18154.981253343561249,
		0.065682337918207449113
	    };
	private final static double[] b = {
		47.20258190468824187,
		976.09855173777669322,
		10260.932208618978205,
		45507.789335026729956
	    };
	
	    private final static double[] c = {
		0.39894151208813466764,
		8.8831497943883759412,
		93.506656132177855979,
		597.27027639480026226,
		2494.5375852903726711,
		6848.1904505362823326,
		11602.651437647350124,
		9842.7148383839780218,
		1.0765576773720192317e-8
	    };
	    private final static double[] d = {
		22.266688044328115691,
		235.38790178262499861,
		1519.377599407554805,
		6485.558298266760755,
		18615.571640885098091,
		34900.952721145977266,
		38912.003286093271411,
		19685.429676859990727
	    };
	    private final static double[] p = {
		0.21589853405795699,
		0.1274011611602473639,
		0.022235277870649807,
		0.001421619193227893466,
		2.9112874951168792e-5,
		0.02307344176494017303
	    };
	    private static final double[] q = {
		1.28426009614491121,
		0.468238212480865118,
		0.0659881378689285515,
		0.00378239633202758244,
		7.29751555083966205e-5
	    };
	    
	private static void pnorm_both(double x, Wrapper<Double> cum, Wrapper<Double> ccum, boolean i_tail, boolean log_p) {
	/* i_tail in {0,1,2} means: "lower", "upper", or "both" :
	   if(lower) return  *cum := P[X <= x]
	   if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
	*/
	    double xden, xnum, temp, del, eps, xsq, y;
/*
	#ifdef NO_DENORMS
	    double min = DBL_MIN;
	#endif
*/
	    int i;
	    boolean lower, upper;

	    /* Consider changing these : */
	    eps = Rmath.DBL_EPSILON * 0.5;

	    /* i_tail in {0,1,2} =^= {lower, upper, both} */
	    lower = !i_tail;
	    upper = i_tail;

	    y = Math.abs(x);
	    if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
		if (y > eps) {
		    xsq = x * x;
		    xnum = a[4] * xsq;
		    xden = xsq;
		    for (i = 0; i < 3; ++i) {
			xnum = (xnum + a[i]) * xsq;
			xden = (xden + b[i]) * xsq;
		    }
		} else xnum = xden = 0.0;

		temp = x * (xnum + a[3]) / (xden + b[3]);
		if(lower)  cum.obj = 0.5 + temp;
		if(upper) ccum.obj = 0.5 - temp;
		if(log_p) {
		    if(lower)  cum.obj = Math.log(cum.obj);
		    if(upper) ccum.obj = Math.log(ccum.obj);
		}
	    }
	    else if (y <= Rmath.M_SQRT_32) {

		/* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

		xnum = c[8] * y;
		xden = y;
		for (i = 0; i < 7; ++i) {
		    xnum = (xnum + c[i]) * y;
		    xden = (xden + d[i]) * y;
		}
		temp = (xnum + c[7]) / (xden + d[7]);		

		xsq = Trunc.trunc(y * SIXTEN) / SIXTEN;
		del = (y - xsq) * (y + xsq);
		if(log_p) {
		    cum.obj = (-xsq * xsq * 0.5) + (-del * 0.5) + Math.log(temp);
		    if((lower && x > 0.) || (upper && x <= 0.))
			  ccum.obj = Log1P.log1p(-Math.exp(-xsq * xsq * 0.5) *
					Math.exp(-del * 0.5) * temp);
		}
		else {
		    cum.obj = Math.exp(-xsq * xsq * 0.5) * Math.exp(-del * 0.5) * temp;
		    ccum.obj = 1.0 - cum.obj;
		}
		if (x > 0.) {/* swap  ccum <--> cum */
		    temp = cum.obj; if(lower) cum.obj = ccum.obj; ccum.obj = temp;
		}
	    }

	/* else	  |x| > sqrt(32) = 5.657 :
	 * the next two case differentiations were really for lower=T, log=F
	 * Particularly	 *not*	for  log_p !

	 * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
	 *
	 * Note that we do want symmetry(0), lower/upper -> hence use y
	 */
	    else if((log_p && y < 1e170) /* avoid underflow below */
		/*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
		 * Then, make use of  Abramowitz & Stegun, 26.2.13, something like

		 xsq = x*x;

		 if(xsq * DBL_EPSILON < 1.)
		    del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
		 else
		    del = 0.;
		 *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + log1p(-del);
		 *ccum = log1p(-exp(*cum)); /.* ~ log(1) = 0 *./

	 	 swap_tail;

		 [Yes, but xsq might be infinite.]

		*/
		    || (lower && -37.5193 < x  &&  x < 8.2924)
		    || (upper && -8.2924  < x  &&  x < 37.5193)
		) {

		/* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
		xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
		xnum = p[5] * xsq;
		xden = xsq;
		for (i = 0; i < 4; ++i) {
		    xnum = (xnum + p[i]) * xsq;
		    xden = (xden + q[i]) * xsq;
		}
		temp = xsq * (xnum + p[4]) / (xden + q[4]);
		temp = (Rmath.M_1_SQRT_2PI - temp) / y;

		xsq = Trunc.trunc(x * SIXTEN) / SIXTEN;
		del = (x - xsq) * (x + xsq);
		if(log_p) {
		    cum.obj = (-xsq * xsq * 0.5) + (-del * 0.5) + Math.log(temp);
		    if((lower && x > 0.) || (upper && x <= 0.))
			  ccum.obj = Log1P.log1p(-Math.exp(-xsq * xsq * 0.5) *
					Math.exp(-del * 0.5) * temp);
		}
		else {
		    cum.obj = Math.exp(-xsq * xsq * 0.5) * Math.exp(-del * 0.5) * temp;
		    ccum.obj = 1.0 - cum.obj;
		}
		if (x > 0.) {/* swap  ccum <--> cum */
		    temp = cum.obj; if(lower) cum.obj = ccum.obj; ccum.obj = temp;
		}
	    } else { /* large x such that probs are 0 or 1 */
		if(x > 0) {	cum.obj = DPQ.RD1(log_p); ccum.obj = DPQ.RD0(log_p);	}
		else {	        cum.obj = DPQ.RD0(log_p); ccum.obj = DPQ.RD1(log_p);	}
	    }

/*
	#ifdef NO_DENORMS
	    /* do not return "denormalized" -- we do in R *  //
	    if(log_p) {
		if(*cum > -min)	 *cum = -0.;
		if(*ccum > -min)*ccum = -0.;
	    }
	    else {
		if(*cum < min)	 *cum = 0.;
		if(*ccum < min)	*ccum = 0.;
	    }
	#endif
*/
	}

}
