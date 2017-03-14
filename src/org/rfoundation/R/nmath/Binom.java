package org.rfoundation.R.nmath;

import org.rfoundation.R.Rmath;

public class Binom {
	/*
	 *  AUTHOR
	 *	Catherine Loader, catherine@research.bell-labs.com.
	 *	October 23, 2000.
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
	 *  DESCRIPTION
	 *	Evaluates the "deviance part"
	 *	bd0(x,M) :=  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =
	 *		  =  x * log(x/M) + M - x
	 *	where M = E[X] = n*p (or = lambda), for	  x, M > 0
	 *
	 *	in a manner that should be stable (with small relative error)
	 *	for all x and M=np. In particular for x/np close to 1, direct
	 *	evaluation fails, and evaluation is based on the Taylor series
	 *	of log((1+v)/(1-v)) with v = (x-np)/(x+np).
	 */
	public static double bd0(double x, double np) {
	    double ej, s, s1, v;
	    int j;

	    if(Double.isInfinite(x) || Double.isInfinite(np) || np == 0.0) {
	    	return Double.NaN;
	    }

	    if (Math.abs(x-np) < 0.1*(x+np)) {
			v = (x-np)/(x+np);
			s = (x-np)*v;/* s using v -- change by MM */
			ej = 2*x*v;
			v = v*v;
			for (j=1; ; j++) { /* Taylor series */
			    ej *= v;
			    s1 = s+ej/((j<<1)+1);
			    if (s1==s) /* last term was effectively 0 */
			    	return(s1);
			    s = s1;
			}
	    }
	    
	    /* else:  | x - np |  is not too small */
	    return(x*Math.log(x/np)+np-x);
	}

	/*
	 * AUTHOR
	 *   Catherine Loader, catherine@research.bell-labs.com.
	 *   October 23, 2000.
	 *
	 *  Merge in to R and further tweaks :
	 *	Copyright (C) 2000, The R Core Development Team
	 *	Copyright (C) 2008, The R Foundation
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
	 *   To compute the binomial probability, call dbinom(x,n,p).
	 *   This checks for argument validity, and calls dbinom_raw().
	 *
	 *   dbinom_raw() does the actual computation; note this is called by
	 *   other functions in addition to dbinom().
	 *     (1) dbinom_raw() has both p and q arguments, when one may be represented
	 *         more accurately than the other (in particular, in df()).
	 *     (2) dbinom_raw() does NOT check that inputs x and n are integers. This
	 *         should be done in the calling function, where necessary.
	 *         -- but is not the case at all when called e.g., from df() or dbeta() !
	 *     (3) Also does not check for 0 <= p <= 1 and 0 <= q <= 1 or NaN's.
	 *         Do this in the calling function.
	 */
	public static double dbinom_raw(int x, int n, double p, double q, boolean give_log)
	{
	    double lf, lc;

	    if (p == 0) {
	    	return((x == 0) ? DPQ.RD1(give_log) : DPQ.RD0(give_log));
	    }
	    
	    if (q == 0) {
	    	return((x == n) ? DPQ.RD1(give_log) : DPQ.RD0(give_log));
	    }

	    if (x == 0) {
			if(n == 0)  {
				return DPQ.RD1(give_log);
			}
			
			lc = (p < 0.1) ? -bd0(n,n*q) - n*p : n* Math.log(q);
			
			return DPQ.R_D_exp(lc, give_log);
	    }
	    
	    if (x == n) {
			lc = (q < 0.1) ? -bd0(n,n*p) - n*q : n * Math.log(p);
			return DPQ.R_D_exp(lc, give_log);
	    }
	    
	    if (x < 0 || x > n) {
	    	return DPQ.RD0(give_log);
	    }

	    /* n*p or n*q can underflow to zero if n and p or q are small.  This
	       used to occur in dbeta, and gives NaN as from R 2.3.0.  */
	    lc = Stirlerr.stirlerr(n) - Stirlerr.stirlerr(x) - 
	    		Stirlerr.stirlerr(n-x) - bd0(x,n*p) - bd0(n-x,n*q);

	    /* f = (M_2PI*x*(n-x))/n; could overflow or underflow */
	    /* Upto R 2.7.1:
	     * lf = log(M_2PI) + log(x) + log(n-x) - log(n);
	     * -- following is much better for  x << n : */
	    lf = Math.log(Rmath.M_2PI) + Math.log(x) + Log1P.log1p(- 1.0 * x/n);

	    return DPQ.R_D_exp(lc - 0.5*lf, give_log);
	}

	public static double dbinom(int x, int n, double p, boolean give_log) {
	    if (p < 0 || p > 1 || n < 0) {
			return Double.NaN;
	    }
	    
	    if (x < 0 || Double.isInfinite(x)) {
	    	return DPQ.RD0(give_log);
	    }
	    
	    return dbinom_raw(x, n, p, 1-p, give_log);
	}
}
