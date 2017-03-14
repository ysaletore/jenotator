package org.rfoundation.R.nmath;

import org.rfoundation.R.Rmath;

public class Gamma {
	/*
	 *  Mathlib : A C Library of Special Functions
	 *  Copyright (C) 1998 Ross Ihaka
	 *  Copyright (C) 2000-2001 The R Development Core Team
	 *  Copyright (C) 2002-2004 The R Foundation
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
	 *    #include <Rmath.h>
	 *    double gammafn(double x);
	 *
	 *  DESCRIPTION
	 *
	 *    This function computes the value of the gamma function.
	 *
	 *  NOTES
	 *
	 *    This function is a translation into C of a Fortran subroutine
	 *    by W. Fullerton of Los Alamos Scientific Laboratory.
	 *    (e.g. http://www.netlib.org/slatec/fnlib/gamma.f)
	 *
	 *    The accuracy of this routine compares (very) favourably
	 *    with those of the Sun Microsystems portable mathematical
	 *    library.
	 *
	 *    MM specialized the case of  n!  for n < 50 - for even better precision
	 */

	private static double[] gamcs = {
		+.8571195590989331421920062399942e-2,
		+.4415381324841006757191315771652e-2,
		+.5685043681599363378632664588789e-1,
		-.4219835396418560501012500186624e-2,
		+.1326808181212460220584006796352e-2,
		-.1893024529798880432523947023886e-3,
		+.3606925327441245256578082217225e-4,
		-.6056761904460864218485548290365e-5,
		+.1055829546302283344731823509093e-5,
		-.1811967365542384048291855891166e-6,
		+.3117724964715322277790254593169e-7,
		-.5354219639019687140874081024347e-8,
		+.9193275519859588946887786825940e-9,
		-.1577941280288339761767423273953e-9,
		+.2707980622934954543266540433089e-10,
		-.4646818653825730144081661058933e-11,
		+.7973350192007419656460767175359e-12,
		-.1368078209830916025799499172309e-12,
		+.2347319486563800657233471771688e-13,
		-.4027432614949066932766570534699e-14,
		+.6910051747372100912138336975257e-15,
		-.1185584500221992907052387126192e-15,
		+.2034148542496373955201026051932e-16,
		-.3490054341717405849274012949108e-17,
		+.5987993856485305567135051066026e-18,
		-.1027378057872228074490069778431e-18,
		+.1762702816060529824942759660748e-19,
		-.3024320653735306260958772112042e-20,
		+.5188914660218397839717833550506e-21,
		-.8902770842456576692449251601066e-22,
		+.1527474068493342602274596891306e-22,
		-.2620731256187362900257328332799e-23,
		+.4496464047830538670331046570666e-24,
		-.7714712731336877911703901525333e-25,
		+.1323635453126044036486572714666e-25,
		-.2270999412942928816702313813333e-26,
		+.3896418998003991449320816639999e-27,
		-.6685198115125953327792127999999e-28,
		+.1146998663140024384347613866666e-28,
		-.1967938586345134677295103999999e-29,
		+.3376448816585338090334890666666e-30,
		-.5793070335782135784625493333333e-31
	    };

	/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
	 * (xmin, xmax) are non-trivial, see ./gammalims.c
	 * xsml = exp(.01)*DBL_MIN
	 * dxrel = sqrt(DBL_EPSILON) = 2 ^ -26
	*/
	private static int ngam = 22;
	private static double xmin = -170.5674972726612;
	private static double xmax = 171.61447887182298;
	private static double xsml = 2.2474362225598545e-308;
	private static double dxrel = 1.490116119384765696e-8;
	
	public static double gammafn(double x) {
		int i, n;
	    double y;
	    double sinpiy, value;
	    
	    if(Double.isNaN(x)) return x;

	    /* If the argument is exactly zero or a negative integer
	     * then return NaN. */
	    if (x == 0 || (x < 0 && x == (long)x)) {
	    	//ML_ERROR(ME_DOMAIN, "gammafn");
	    	return Double.NaN;
	    }

	    y = Math.abs(x);

	    if (y <= 10) {

		/* Compute gamma(x) for -10 <= x <= 10
		 * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
		 * first of all. */

		n = (int) x;
		if(x < 0) --n;
		y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
		--n;
		value = Chebyshev.chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
		if (n == 0)
		    return value;/* x = 1.dddd = 1+y */

		if (n < 0) {
		    /* compute gamma(x) for -10 <= x < 1 */

		    /* exact 0 or "-n" checked already above */

		    /* The answer is less than half precision */
		    /* because x too near a negative integer. */
		    if (x < -0.5 && Math.abs(x - (int)(x - 0.5) / x) < dxrel) {
		    	// TODO ML_ERROR(ME_PRECISION, "gammafn");
		    }

		    /* The argument is so close to 0 that the result would overflow. */
		    if (y < xsml) {
		    	//	TODO ML_ERROR(ME_RANGE, "gammafn");
				if(x > 0) return Double.POSITIVE_INFINITY;
				else return Double.NEGATIVE_INFINITY;
		    }

		    n = -n;

		    for (i = 0; i < n; i++) {
			value /= (x + i);
		    }
		    return value;
		}
		else {
		    /* gamma(x) for 2 <= x <= 10 */

		    for (i = 1; i <= n; i++) {
			value *= (y + i);
		    }
		    return value;
		}
	    }
	    else {
		/* gamma(x) for	 y = |x| > 10. */

		if (x > xmax) {			/* Overflow */
		    // TODO ML_ERROR(ME_RANGE, "gammafn");
		    return Double.POSITIVE_INFINITY;
		}

		if (x < xmin) {			/* Underflow */
		    // TODO ML_ERROR(ME_UNDERFLOW, "gammafn");
		    return 0.;
		}

		if(y <= 50 && y == (int)y) { /* compute (n - 1)! */
		    value = 1.;
		    for (i = 2; i < y; i++) value *= i;
		}
		else { /* normal case */
		    value = Math.exp((y - 0.5) * Math.log(y) - y + Rmath.M_LN_SQRT_2PI +
				((2*y == (int)2*y)? Stirlerr.stirlerr(y) : LGammaCor.lgammacor(y)));
		}
		if (x > 0)
		    return value;

		if (Math.abs((x - (int)(x - 0.5))/x) < dxrel){

		    /* The answer is less than half precision because */
		    /* the argument is too near a negative integer. */

		    // FIXME ML_ERROR(ME_PRECISION, "gammafn");
		}

		sinpiy = Math.sin(Rmath.M_PI * y);
		if (sinpiy == 0) {		/* Negative integer arg - overflow */
		    // FIXME ML_ERROR(ME_RANGE, "gammafn");
		    return Double.POSITIVE_INFINITY;
		}

		return -Rmath.M_PI / (y * sinpiy * value);
	    }
	}

	/*
	 *  AUTHOR
	 *    Catherine Loader, catherine@research.bell-labs.com.
	 *    October 23, 2000.
	 *
	 *  Merge in to R:
	 *	Copyright (C) 2000 The R Core Development Team
	 *	Copyright (C) 2004 The R Foundation
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
	 *
	 * DESCRIPTION
	 *
	 *   Computes the density of the gamma distribution,
	 *
	 *                   1/s (x/s)^{a-1} exp(-x/s)
	 *        p(x;a,s) = -----------------------
	 *                            (a-1)!
	 *
	 *   where `s' is the scale (= 1/lambda in other parametrizations)
	 *     and `a' is the shape parameter ( = alpha in other contexts).
	 *
	 * The old (R 1.1.1) version of the code is available via `#define D_non_pois'
	 */
	public static double dgamma(double x, double shape, double scale, boolean give_log) {
	    double pr;

	    if (shape < 0 || scale <= 0)  {
	    	return Double.NaN;
	    }
	    
	    if (x < 0) {
	    	return DPQ.RD0(give_log);
	    }
	    
	    if (shape == 0) { /* point mass at 0 */
	    	return (x == 0) ? Double.POSITIVE_INFINITY : DPQ.RD0(give_log);
	    }
	    
	    if (x == 0) {
			if (shape < 1) {
				return Double.POSITIVE_INFINITY;
			}
			
			if (shape > 1) {
				return DPQ.RD0(give_log);
			}
			
			/* else */
			return give_log ? -Math.log(scale) : 1 / scale;
	    }

	    if (shape < 1) {
	    	pr = dpois_raw(shape, x/scale, give_log);
	    	return give_log ?  pr + Math.log(shape/x) : pr*shape/x;
	    }
	    /* else  shape >= 1 */
	    pr = dpois_raw(shape-1, x/scale, give_log);
	    return give_log ? pr - Math.log(scale) : pr/scale;
	}
	
	/*
	 *  AUTHOR
	 *    Catherine Loader, catherine@research.bell-labs.com.
	 *    October 23, 2000.
	 *
	 *  Merge in to R:
	 *	Copyright (C) 2000, The R Core Development Team
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
	 *
	 * DESCRIPTION
	 *
	 *    dpois() checks argument validity and calls dpois_raw().
	 *
	 *    dpois_raw() computes the Poisson probability  lb^x exp(-lb) / x!.
	 *      This does not check that x is an integer, since dgamma() may
	 *      call this with a fractional x argument. Any necessary argument
	 *      checks should be done in the calling function.
	 *
	 */

	private static double dpois_raw(double x, double lambda, boolean give_log) {
	    /*       x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
	        lambda >= 0
	    */
	    if (lambda == 0) return( (x == 0) ? DPQ.RD1(give_log) : DPQ.RD0(give_log) );
	    if (Double.isInfinite(lambda)) return DPQ.RD0(give_log);
	    if (x < 0) return( DPQ.RD0(give_log) );
	    if (x <= lambda * Double.MIN_VALUE) return(DPQ.R_D_exp(-lambda, give_log) );
	    if (lambda < x * Double.MIN_VALUE) return DPQ.R_D_exp(-lambda + x*Math.log(lambda) -LGamma.lgammafn(x+1), give_log);
	    return DPQ.R_D_fexp( Rmath.M_2PI*x, -Stirlerr.stirlerr(x)-Binom.bd0(x,lambda), give_log);
	}

	private static double dpois(double x, double lambda, boolean give_log) {
	    if (lambda < 0) {
	    	return Double.NaN;
	    }
	    
	    // R_D_nonint_check(x);
	    if (x < 0 || Double.isInfinite(x))
		return DPQ.RD0(give_log);

	    x = Math.floor(x+0.5);

	    return( dpois_raw(x,lambda,give_log) );
	}

	/*
	 *  Mathlib : A C Library of Special Functions
	 *  Copyright (C) 2005-6 Morten Welinder <terra@gnome.org>
	 *  Copyright (C) 2005-10 The R Foundation
	 *  Copyright (C) 2006-10 The R Core Development Team
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
	 *	#include <Rmath.h>
	 *
	 *	double pgamma (double x, double alph, double scale,
	 *		       boolean lower_tail, boolean log_p)
	 *
	 *	double log1pmx	(double x)
	 *	double lgamma1p (double a)
	 *
	 *	double logspace_add (double logx, double logy)
	 *	double logspace_sub (double logx, double logy)
	 *
	 *
	 *  DESCRIPTION
	 *
	 *	This function computes the distribution function for the
	 *	gamma distribution with shape parameter alph and scale parameter
	 *	scale.	This is also known as the incomplete gamma function.
	 *	See Abramowitz and Stegun (6.5.1) for example.
	 *
	 *  NOTES
	 *
	 *	Complete redesign by Morten Welinder, originally for Gnumeric.
	 *	Improvements (e.g. "while NEEDED_SCALE") by Martin Maechler
	 *	The old version can be activated by compiling with -DR_USE_OLD_PGAMMA
	 *
	 *  REFERENCES
	 *
	 */

	/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
	public static double SQR(double x ) {
		return x * x;
	}
	
	private static final double scalefactor = SQR(SQR(SQR(4294967296.0)));

	/* If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ]	 =~=  -x */
	private static final double M_cutoff = Rmath.M_LN2 * Double.MAX_EXPONENT / Rmath.DBL_EPSILON;/*=3.196577e18*/

	/* Continued fraction for calculation of
	 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
	 *
	 * auxilary in log1pmx() and lgamma1p()
	 */
	private static double
	logcf (double x, double i, double d,
	       double eps /* ~ relative tolerance */)
	{
	    double c1 = 2 * d;
	    double c2 = i + d;
	    double c4 = c2 + d;
	    double a1 = c2;
	    double b1 = i * (c2 - i * x);
	    double b2 = d * d * x;
	    double a2 = c4 * c2 - b2;

	    b2 = c4 * b1 - i * b2;

	    while (Math.abs(a2 * b1 - a1 * b2) > Math.abs(eps * b1 * b2)) {
		double c3 = c2*c2*x;
		c2 += d;
		c4 += d;
		a1 = c4 * a2 - c3 * a1;
		b1 = c4 * b2 - c3 * b1;

		c3 = c1 * c1 * x;
		c1 += d;
		c4 += d;
		a2 = c4 * a1 - c3 * a2;
		b2 = c4 * b1 - c3 * b2;

		if (Math.abs (b2) > scalefactor) {
		    a1 /= scalefactor;
		    b1 /= scalefactor;
		    a2 /= scalefactor;
		    b2 /= scalefactor;
		} else if (Math.abs (b2) < 1 / scalefactor) {
		    a1 *= scalefactor;
		    b1 *= scalefactor;
		    a2 *= scalefactor;
		    b2 *= scalefactor;
		}
	    }

	    return a2 / b2;
	}

	private static final double minLog1Value = -0.79149064;
	private static final double two = 2;
	private static final double tol_logcf = 1e-14;
	/* Accurate calculation of log(1+x)-x, particularly for small x.  */
	private static double log1pmx (double x) {
	    if (x > 1 || x < minLog1Value)
		return Log1P.log1p(x) - x;
	    else { /* -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
		    * log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
		    * ---------------------------------------------
		    * S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
		   */
		double r = x / (2 + x), y = r * r;
		if (Math.abs(x) < 1e-2) {
		    
		    return r * ((((two / 9 * y + two / 7) * y + two / 5) * y +
				    two / 3) * y - x);
		} else {
		    
		    return r * (2 * y * logcf (y, 3, 2, tol_logcf) - x);
		}
	    }
	}


	private static final double eulers_const =	 0.5772156649015328606065120900824024;
	 /* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
    private static final int N = 40;
    private static final double[] coeffs = {
	0.3224670334241132182362075833230126e-0,/* = (zeta(2)-1)/2 */
	0.6735230105319809513324605383715000e-1,/* = (zeta(3)-1)/3 */
	0.2058080842778454787900092413529198e-1,
	0.7385551028673985266273097291406834e-2,
	0.2890510330741523285752988298486755e-2,
	0.1192753911703260977113935692828109e-2,
	0.5096695247430424223356548135815582e-3,
	0.2231547584535793797614188036013401e-3,
	0.9945751278180853371459589003190170e-4,
	0.4492623673813314170020750240635786e-4,
	0.2050721277567069155316650397830591e-4,
	0.9439488275268395903987425104415055e-5,
	0.4374866789907487804181793223952411e-5,
	0.2039215753801366236781900709670839e-5,
	0.9551412130407419832857179772951265e-6,
	0.4492469198764566043294290331193655e-6,
	0.2120718480555466586923135901077628e-6,
	0.1004322482396809960872083050053344e-6,
	0.4769810169363980565760193417246730e-7,
	0.2271109460894316491031998116062124e-7,
	0.1083865921489695409107491757968159e-7,
	0.5183475041970046655121248647057669e-8,
	0.2483674543802478317185008663991718e-8,
	0.1192140140586091207442548202774640e-8,
	0.5731367241678862013330194857961011e-9,
	0.2759522885124233145178149692816341e-9,
	0.1330476437424448948149715720858008e-9,
	0.6422964563838100022082448087644648e-10,
	0.3104424774732227276239215783404066e-10,
	0.1502138408075414217093301048780668e-10,
	0.7275974480239079662504549924814047e-11,
	0.3527742476575915083615072228655483e-11,
	0.1711991790559617908601084114443031e-11,
	0.8315385841420284819798357793954418e-12,
	0.4042200525289440065536008957032895e-12,
	0.1966475631096616490411045679010286e-12,
	0.9573630387838555763782200936508615e-13,
	0.4664076026428374224576492565974577e-13,
	0.2273736960065972320633279596737272e-13,
	0.1109139947083452201658320007192334e-13/* = (zeta(40+1)-1)/(40+1) */
    };

    private static final double c = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
    // private static final double tol_logcf = 1e-14;
	/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
	private static double lgamma1p (double a) {
	    double lgam;
	    int i;

	    if (Math.abs (a) >= 0.5)
		return LGamma.lgammafn (a + 1);

	    /* Abramowitz & Stegun 6.1.33 : for |x| < 2,
	     * <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
	     * where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
	     *
	     * Here, another convergence acceleration trick is used to compute
	     * lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n
	     */
	    lgam = c * logcf(-a / 2, N + 2, 1, tol_logcf);
	    for (i = N - 1; i >= 0; i--)
		lgam = coeffs[i] - a * lgam;

	    return (a * lgam - eulers_const) * a - log1pmx (a);
	} /* lgamma1p */



	/*
	 * Compute the log of a sum from logs of terms, i.e.,
	 *
	 *     log (exp (logx) + exp (logy))
	 *
	 * without causing overflows and without throwing away large handfuls
	 * of accuracy.
	 */
	private static double logspace_add (double logx, double logy)
	{
	    return Math.max (logx, logy) + Log1P.log1p (Math.exp (-Math.abs (logx - logy)));
	}


	/*
	 * Compute the log of a difference from logs of terms, i.e.,
	 *
	 *     log (exp (logx) - exp (logy))
	 *
	 * without causing overflows and without throwing away large handfuls
	 * of accuracy.
	 */
	double logspace_sub (double logx, double logy)
	{
	    return logx + Log1P.log1p (-Math.exp (logy - logx));
	}


	// #ifndef R_USE_OLD_PGAMMA

	/* dpois_wrap (x_P_1,  lambda, g_log) ==
	 *   dpois (x_P_1 - 1, lambda, g_log) :=  exp(-L)  L^k / gamma(k+1) ,  k := x_P_1 - 1
	*/
	private static double
	dpois_wrap (double x_plus_1, double lambda, boolean give_log) {
	    if (Double.isInfinite(lambda))
		return DPQ.RD0(give_log);
	    if (x_plus_1 > 1)
		return dpois_raw (x_plus_1 - 1, lambda, give_log);
	    if (lambda > Math.abs(x_plus_1 - 1) * M_cutoff)
		return DPQ.R_D_exp(-lambda - LGamma.lgammafn(x_plus_1), give_log);
	    else {
		double d = dpois_raw (x_plus_1, lambda, give_log);

		return give_log
		    ? d + Math.log (x_plus_1 / lambda)
		    : d * (x_plus_1 / lambda);
	    }
	}

	/*
	 * Abramowitz and Stegun 6.5.29 [right]
	 */
	static double
	pgamma_smallx (double x, double alph, boolean lower_tail, boolean log_p)
	{
	    double sum = 0, c = alph, n = 0, term;

	    /*
	     * Relative to 6.5.29 all terms have been multiplied by alph
	     * and the first, thus being 1, is omitted.
	     */

	    do {
		n++;
		c *= -x / n;
		term = c / (alph + n);
		sum += term;
	    } while (Math.abs (term) > Rmath.DBL_EPSILON * Math.abs (sum));
	    
	    if (lower_tail) {
		double f1 = log_p ? Log1P.log1p (sum) : 1 + sum;
		double f2;
		if (alph > 1) {
		    f2 = dpois_raw (alph, x, log_p);
		    f2 = log_p ? f2 + x : f2 * Math.exp (x);
		} else if (log_p)
		    f2 = alph * Math.log (x) - lgamma1p (alph);
		else
		    f2 = Math.pow (x, alph) / Math.exp (lgamma1p (alph));
		return log_p ? f1 + f2 : f1 * f2;
	    } else {
		double lf2 = alph * Math.log (x) - lgamma1p (alph);
		if (log_p)
		    return DPQ.R_Log1_Exp (Log1P.log1p (sum) + lf2);
		else {
		    double f1m1 = sum;
		    double f2m1 = ExpM1.expm1 (lf2);
		    return -(f1m1 + f2m1 + f1m1 * f2m1);
		}
	    }
	} /* pgamma_smallx() */

	private static double pd_upper_series (double x, double y, boolean log_p) {
	    double term = x / y;
	    double sum = term;

	    do {
		y++;
		term *= x / y;
		sum += term;
	    } while (term > sum * Rmath.DBL_EPSILON);

	    /* sum =  \sum_{n=1}^ oo  x^n     / (y*(y+1)*...*(y+n-1))
	     *	   =  \sum_{n=0}^ oo  x^(n+1) / (y*(y+1)*...*(y+n))
	     *	   =  x/y * (1 + \sum_{n=1}^oo	x^n / ((y+1)*...*(y+n)))
	     *	   ~  x/y +  o(x/y)   {which happens when alph -> Inf}
	     */
	    return log_p ? Math.log (sum) : sum;
	}

	/* Continued fraction for calculation of
	 *    scaled upper-tail F_{gamma}
	 *  ~=  (y / d) * [1 +  (1-y)/d +  O( ((1-y)/d)^2 ) ]
	 */
	private static int max_it = 200000;
	private static double pd_lower_cf (double y, double d)
	{
	    double f= 0.0 /* -Wall */, of, f0;
	    double i, c2, c3, c4,  a1, b1,  a2, b2;

	    if (y == 0) return 0;

	    f0 = y/d;
	    /* Needed, e.g. for  pgamma(10^c(100,295), shape= 1.1, log=TRUE): */
	    if(Math.abs(y - 1) < Math.abs(d) * Rmath.DBL_EPSILON) { /* includes y < d = Inf */
		return (f0);
	    }

	    if(f0 > 1.) f0 = 1.;
	    c2 = y;
	    c4 = d; /* original (y,d), *not* potentially scaled ones!*/

	    a1 = 0; b1 = 1;
	    a2 = y; b2 = d;

	    while (b2 > scalefactor) {
		    a1 /= scalefactor;
		    b1 /= scalefactor;
		    a2 /= scalefactor;
		    b2 /= scalefactor;
		}

	    i = 0; of = -1.; /* far away */
	    while (i < max_it) {

		i++;	c2--;	c3 = i * c2;	c4 += 2;
		/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd */
		a1 = c4 * a2 + c3 * a1;
		b1 = c4 * b2 + c3 * b1;

		i++;	c2--;	c3 = i * c2;	c4 += 2;
		/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even */
		a2 = c4 * a1 + c3 * a2;
		b2 = c4 * b1 + c3 * b2;

		if (b2 > scalefactor) {
		    a1 /= scalefactor;
		    b1 /= scalefactor;
		    a2 /= scalefactor;
		    b2 /= scalefactor;
		}

		if (b2 != 0) {
		    f = a2 / b2;
	 	    /* convergence check: relative; "absolute" for very small f : */
		    if (Math.abs (f - of) <= Rmath.DBL_EPSILON * Math.max(f0, Math.abs(f))) {
			return f;
		    }
		    of = f;
		}
	    }

	  // FIXME  MATHLIB_WARNING(" ** NON-convergence in pgamma()'s pd_lower_cf() f= %g.\n",
	  //	    f);
	    return f;/* should not happen ... */
	} /* pd_lower_cf() */

	static double
	pd_lower_series (double lambda, double y)
	{
	    double term = 1, sum = 0;

	    while (y >= 1 && term > sum * Rmath.DBL_EPSILON) {
		term *= y / lambda;
		sum += term;
		y--;
	    }
	    /* sum =  \sum_{n=0}^ oo  y*(y-1)*...*(y - n) / lambda^(n+1)
	     *	   =  y/lambda * (1 + \sum_{n=1}^Inf  (y-1)*...*(y-n) / lambda^n)
	     *	   ~  y/lambda + o(y/lambda)
	     */

	    if (y != Math.floor (y)) {
		/*
		 * The series does not converge as the terms start getting
		 * bigger (besides flipping sign) for y < -lambda.
		 */
		double f;

		/* FIXME: in quite few cases, adding  term*f  has no effect (f too small)
		 *	  and is unnecessary e.g. for pgamma(4e12, 121.1) */
		f = pd_lower_cf (y, lambda + 1 - y);
		
		sum += term * f;
	    }

	    return sum;
	} /* pd_lower_series() */

	/*
	 * Compute the following ratio with higher accuracy that would be had
	 * from doing it directly.
	 *
	 *		 dnorm (x, 0, 1, FALSE)
	 *	   ----------------------------------
	 *	   pnorm (x, 0, 1, lower_tail, FALSE)
	 *
	 * Abramowitz & Stegun 26.2.12
	 */
	private static double dpnorm (double x, boolean lower_tail, double lp)
	{
	    /*
	     * So as not to repeat a pnorm call, we expect
	     *
	     *	 lp == pnorm (x, 0, 1, lower_tail, TRUE)
	     *
	     * but use it only in the non-critical case where either x is small
	     * or p==exp(lp) is close to 1.
	     */

	    if (x < 0) {
		x = -x;
		lower_tail = !lower_tail;
	    }

	    if (x > 10 && !lower_tail) {
		double term = 1 / x;
		double sum = term;
		double x2 = x * x;
		double i = 1;

		do {
		    term *= -i / x2;
		    sum += term;
		    i += 2;
		} while (Math.abs (term) > Rmath.DBL_EPSILON * sum);

		return 1 / sum;
	    } else {
		double d = Norm.dnorm (x, 0, 1, false);
		return d / Math.exp (lp);
	    }
	}

	/*
	 * Asymptotic expansion to calculate the probability that Poisson variate
	 * has value <= x.
	 * Various assertions about this are made (without proof) at
	 * http://members.aol.com/iandjmsmith/PoissonApprox.htm
	 */
	private static final double[] coefs_a = {
		-1e99, /* placeholder used for 1-indexing */
		2/3.,
		-4/135.,
		8/2835.,
		16/8505.,
		-8992/12629925.,
		-334144/492567075.,
		698752/1477701225.
	    };
	
	private static final double[] coefs_b = {
		-1e99, /* placeholder */
		1/12.,
		1/288.,
		-139/51840.,
		-571/2488320.,
		163879/209018880.,
		5246819/75246796800.,
		-534703531/902961561600.
	    };
	private static double ppois_asymp (double x, double lambda, boolean lower_tail, boolean log_p) {
	    double elfb, elfb_term;
	    double res12, res1_term, res1_ig, res2_term, res2_ig;
	    double dfm, pt_, s2pt, f, np;
	    int i;

	    dfm = lambda - x;
	    /* If lambda is large, the distribution is highly concentrated
	       about lambda.  So representation error in x or lambda can lead
	       to arbitrarily large values of pt_ and hence divergence of the
	       coefficients of this approximation.
	    */
	    pt_ = - log1pmx (dfm / x);
	    s2pt = Math.sqrt (2 * x * pt_);
	    if (dfm < 0) s2pt = -s2pt;

	    res12 = 0;
	    res1_ig = res1_term = Math.sqrt (x);
	    res2_ig = res2_term = s2pt;
	    for (i = 1; i < 8; i++) {
		res12 += res1_ig * coefs_a[i];
		res12 += res2_ig * coefs_b[i];
		res1_term *= pt_ / i ;
		res2_term *= 2 * pt_ / (2 * i + 1);
		res1_ig = res1_ig / x + res1_term;
		res2_ig = res2_ig / x + res2_term;
	    }

	    elfb = x;
	    elfb_term = 1;
	    for (i = 1; i < 8; i++) {
		elfb += elfb_term * coefs_b[i];
		elfb_term /= x;
	    }
	    if (!lower_tail) elfb = -elfb;

	    f = res12 / elfb;

	    np = Norm.pnorm (s2pt, 0.0, 1.0, !lower_tail, log_p);

	    if (log_p) {
		double n_d_over_p = dpnorm (s2pt, !lower_tail, np);
		return np + Log1P.log1p (f * n_d_over_p);
	    } else {
		double nd = Norm.dnorm(s2pt, 0., 1., log_p);

		return np + f * nd;
	    }
	} /* ppois_asymp() */


	private static double pgamma_raw (double x, double alph, boolean lower_tail, boolean log_p) {
	/* Here, assume that  (x,alph) are not NA  &  alph > 0 . */

	    double res;

	    if(x <= 0.) {
	    	return DPQ.R_DT_0(lower_tail, log_p);
	    }
	    
	    if(x >= Double.POSITIVE_INFINITY) {
	    	return DPQ.R_DT_1(lower_tail, log_p);
	    }

	    if (x < 1) {
		res = pgamma_smallx (x, alph, lower_tail, log_p);
	    } else if (x <= alph - 1 && x < 0.8 * (alph + 50)) {
		/* incl. large alph compared to x */
		double sum = pd_upper_series (x, alph, log_p);/* = x/alph + o(x/alph) */
		double d = dpois_wrap (alph, x, log_p);

		if (!lower_tail)
		    res = log_p
			? DPQ.R_Log1_Exp (d + sum)
			: 1 - d * sum;
		else
		    res = log_p ? sum + d : sum * d;
	    } else if (alph - 1 < x && alph < 0.8 * (x + 50)) {
		/* incl. large x compared to alph */
		double sum;
		double d = dpois_wrap (alph, x, log_p);

		if (alph < 1) {
		    if (x * Rmath.DBL_EPSILON > 1 - alph)
			sum = DPQ.RD1(log_p);
		    else {
			double f = pd_lower_cf (alph, x - (alph - 1)) * x / alph;
			/* = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1) */
			sum = log_p ? Math.log (f) : f;
		    }
		} else {
		    sum = pd_lower_series (x, alph - 1);/* = (alph-1)/x + o((alph-1)/x) */
		    sum = log_p ? Log1P.log1p (sum) : 1 + sum;
		}

		if (!lower_tail)
		    res = log_p ? sum + d : sum * d;
		else
		    res = log_p
			? DPQ.R_Log1_Exp (d + sum)
			: 1 - d * sum;
	    } else { /* x >= 1 and x fairly near alph. */
		res = ppois_asymp (alph - 1, x, !lower_tail, log_p);
	    }

	    /*
	     * We lose a fair amount of accuracy to underflow in the cases
	     * where the final result is very close to DBL_MIN.	 In those
	     * cases, simply redo via log space.
	     */
	    if (!log_p && res < Double.MIN_VALUE / Rmath.DBL_EPSILON) {
		/* with(.Machine, double.xmin / double.eps) #|-> 1.002084e-292 */
		return Math.exp (pgamma_raw (x, alph, lower_tail, true));
	    } else
		return res;
	}

	public static double pgamma(double x, double alph, double scale, boolean lower_tail, boolean log_p) {
	    if(alph < 0. || scale <= 0.)
			return Double.NaN;
	    x /= scale;

	    if(alph == 0.) /* limit case; useful e.g. in pnchisq() */
	    	return (x <= 0) ? DPQ.R_DT_0(lower_tail, log_p): DPQ.R_DT_1(lower_tail, log_p); /* <= assert  pgamma(0,0) ==> 0 */
	    return pgamma_raw (x, alph, lower_tail, log_p);
	}
	/* From: terra@gnome.org (Morten Welinder)
	 * To: R-bugs@biostat.ku.dk
	 * Cc: maechler@stat.math.ethz.ch
	 * Subject: Re: [Rd] pgamma discontinuity (PR#7307)
	 * Date: Tue, 11 Jan 2005 13:57:26 -0500 (EST)

	 * this version of pgamma appears to be quite good and certainly a vast
	 * improvement over current R code.  (I last looked at 2.0.1)  Apart from
	 * type naming, this is what I have been using for Gnumeric 1.4.1.

	 * This could be included into R as-is, but you might want to benefit from
	 * making logcf, log1pmx, lgamma1p, and possibly logspace_add/logspace_sub
	 * available to other parts of R.

	 * MM: I've not (yet?) taken  logcf(), but the other four
	 */

/*
	#else
	/* R_USE_OLD_PGAMMA */
	/*
	 *  Copyright (C) 1998		Ross Ihaka
	 *  Copyright (C) 1999-2000	The R Development Core Team
	 *  Copyright (C) 2003-2004	The R Foundation
	 *  based on AS 239 (C) 1988 Royal Statistical Society
	 *
	 *  ................
	 *
	 *  NOTES
	 *
	 *	This function is an adaptation of Algorithm 239 from the
	 *	Applied Statistics Series.  The algorithm is faster than
	 *	those by W. Fullerton in the FNLIB library and also the
	 *	TOMS 542 alorithm of W. Gautschi.  It provides comparable
	 *	accuracy to those algorithms and is considerably simpler.
	 *
	 *  REFERENCES
	 *
	 *	Algorithm AS 239, Incomplete Gamma Function
	 *	Applied Statistics 37, 1988.
	 * //

	/* now would need this here: * //
	double attribute_hidden pgamma_raw(x, alph, lower_tail, log_p) {
	    return pgamma(x, alph, 1, lower_tail, log_p);
	}

	double pgamma(double x, double alph, double scale, boolean lower_tail, boolean log_p)
	{
	    final static double
		xbig = 1.0e+8,
		xlarge = 1.0e+37,

		/* normal approx. for alph > alphlimit * //
		alphlimit = 1e5;/* was 1000. till R.1.8.x * //

	    double pn1, pn2, pn3, pn4, pn5, pn6, arg, a, b, c, an, osum, sum;
	    long n;
	    int pearson;

	    /* check that we have valid values for x and alph * //

	#ifdef IEEE_754
	    if (ISNAN(x) || ISNAN(alph) || ISNAN(scale))
		return x + alph + scale;
	#endif
	#ifdef DEBUG_p
	    REprintf("pgamma(x=%4g, alph=%4g, scale=%4g): ",x,alph,scale);
	#endif
	    if(alph <= 0. || scale <= 0.)
		ML_ERR_return_NAN;

	    x /= scale;
	#ifdef DEBUG_p
	    REprintf("-> x=%4g; ",x);
	#endif
	#ifdef IEEE_754
	    if (ISNAN(x)) /* eg. original x = scale = Inf *  //
		return x;
	#endif
	    if (x <= 0.)
		return R_DT_0;

	#define USE_PNORM \
	    pn1 = sqrt(alph) * 3. * (pow(x/alph, 1./3.) + 1. / (9. * alph) - 1.); \
	    return pnorm(pn1, 0., 1., lower_tail, log_p);

	    if (alph > alphlimit) { /* use a normal approximation *  //
		USE_PNORM;
	    }

	    if (x > xbig * alph) {
		if (x > DBL_MAX * alph)
		    /* if x is extremely large __compared to alph__ then return 1 *  //
		    return R_DT_1;
		else { /* this only "helps" when log_p = TRUE *  //
		    USE_PNORM;
		}
	    }

	    if (x <= 1. || x < alph) {

		pearson = 1;/* use pearson's series expansion. *  //

		arg = alph * log(x) - x - lgammafn(alph + 1.);
	#ifdef DEBUG_p
		REprintf("Pearson  arg=%g ", arg);
	#endif
		c = 1.;
		sum = 1.;
		a = alph;
		do {
		    a += 1.;
		    c *= x / a;
		    sum += c;
		} while (c > DBL_EPSILON * sum);
	    }
	    else { /* x >= max( 1, alph) * //

		pearson = 0;/* use a continued fraction expansion *  //

		arg = alph * log(x) - x - lgammafn(alph);
	#ifdef DEBUG_p
		REprintf("Cont.Fract. arg=%g ", arg);
	#endif
		a = 1. - alph;
		b = a + x + 1.;
		pn1 = 1.;
		pn2 = x;
		pn3 = x + 1.;
		pn4 = x * b;
		sum = pn3 / pn4;
		for (n = 1; ; n++) {
		    a += 1.;/* =   n+1 -alph *  //
		    b += 2.;/* = 2(n+1)-alph+x *  //
		    an = a * n;
		    pn5 = b * pn3 - an * pn1;
		    pn6 = b * pn4 - an * pn2;
		    if (Math.abs(pn6) > 0.) {
			osum = sum;
			sum = pn5 / pn6;
			if (Math.abs(osum - sum) <= DBL_EPSILON * fmin2(1., sum))
			    break;
		    }
		    pn1 = pn3;
		    pn2 = pn4;
		    pn3 = pn5;
		    pn4 = pn6;
		    if (Math.abs(pn5) >= xlarge) {
			/* re-scale the terms in continued fraction if they are large *  //
	#ifdef DEBUG_p
			REprintf(" [r] ");
	#endif
			pn1 /= xlarge;
			pn2 /= xlarge;
			pn3 /= xlarge;
			pn4 /= xlarge;
		    }
		}
	    }

	    arg += log(sum);

	    lower_tail = (lower_tail == pearson);

	    if (log_p && lower_tail)
		return(arg);
	    /* else */
	    /* sum = exp(arg); and return   if(lower_tail) sum	else 1-sum : *  //
	    return (lower_tail) ? exp(arg) : (log_p ? R_Log1_Exp(arg) : -expm1(arg));
	}

	#endif
	/* R_USE_OLD_PGAMMA */
}
