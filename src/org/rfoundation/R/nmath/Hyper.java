package org.rfoundation.R.nmath;

import org.rfoundation.R.Rmath;

public class Hyper {
	/*
	 *  AUTHOR
	 *    Catherine Loader, catherine@research.bell-labs.com.
	 *    October 23, 2000.
	 *
	 *  Merge in to R:
	 *	Copyright (C) 2000, 2001 The R Core Development Team
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
	 *    Given a sequence of r successes and b failures, we sample n (\le b+r)
	 *    items without replacement. The hypergeometric probability is the
	 *    probability of x successes:
	 *
	 *		       choose(r, x) * choose(b, n-x)
	 *	p(x; r,b,n) =  -----------------------------  =
	 *			       choose(r+b, n)
	 *
	 *		      dbinom(x,r,p) * dbinom(n-x,b,p)
	 *		    = --------------------------------
	 *			       dbinom(n,r+b,p)
	 *
	 *    for any p. For numerical stability, we take p=n/(r+b); with this choice,
	 *    the denominator is not exponentially small.
	 */

	public static double dhyper(int x, int r, int b, int n, boolean give_log) {
	    double p, q, p1, p2, p3;

	    if(r < 0 || b < 0 || n < 0 || n > r + b) {
	    	return Double.NaN;
	    }
	    
	    if(x < 0) {
	    	return DPQ.RD0(give_log);
	    }

	    if (n < x || r < x || n - x > b) return(DPQ.RD0(give_log));
	    if (n == 0) return((x == 0) ? DPQ.RD1(give_log) : DPQ.RD0(give_log));

	    p = ((double)n)/((double)(r+b));
	    q = ((double)(r+b-n))/((double)(r+b));

	    p1 = Binom.dbinom_raw(x, r, p,q,give_log);
	    p2 = Binom.dbinom_raw(n-x,b, p,q,give_log);
	    p3 = Binom.dbinom_raw(n,r+b, p,q,give_log);

	    return( (give_log) ? p1 + p2 - p3 : p1*p2/p3 );
	}
	
	/*
	 *  Mathlib : A C Library of Special Functions
	 *  Copyright (C) 1998 Ross Ihaka
	 *  Copyright (C) 1999-2000  The R Development Core Team
	 *  Copyright (C) 2004	     Morten Welinder
	 *  Copyright (C) 2004	     The R Foundation
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
	 *  DESCRIPTION
	 *
	 *	The distribution function of the hypergeometric distribution.
	 *
	 * Current implementation based on posting
	 * From: Morten Welinder <terra@gnome.org>
	 * Cc: R-bugs@biostat.ku.dk
	 * Subject: [Rd] phyper accuracy and efficiency (PR#6772)
	 * Date: Thu, 15 Apr 2004 18:06:37 +0200 (CEST)
	 ......

	 The current version has very serious cancellation issues.  For example,
	 if you ask for a small right-tail you are likely to get total cancellation.
	 For example,  phyper(59, 150, 150, 60, FALSE, FALSE) gives 6.372680161e-14.
	 The right answer is dhyper(0, 150, 150, 60, FALSE) which is 5.111204798e-22.

	 phyper is also really slow for large arguments.

	 Therefore, I suggest using the code below. This is a sniplet from Gnumeric ...
	 The code isn't perfect.  In fact, if  x*(NR+NB)  is close to	n*NR,
	 then this code can take a while. Not longer than the old code, though.

	 -- Thanks to Ian Smith for ideas.
	*/

	public static double pdhyper (double x, double NR, double NB, double n, boolean log_p) {
	/*
	 * Calculate
	 *
	 *	    phyper (x, NR, NB, n, TRUE, FALSE)
	 *   [log]  ----------------------------------
	 *	       dhyper (x, NR, NB, n, FALSE)
	 *
	 * without actually calling phyper.  This assumes that
	 *
	 *     x * (NR + NB) <= n * NR
	 *
	 */
	    double sum = 0;
	    double term = 1;

	    while (x > 0 && term >= Rmath.DBL_EPSILON * sum) {
			term *= x * (NB - n + x) / (n + 1 - x) / (NR + 1 - x);
			sum += term;
			x--;
	    }

	    return log_p ? Log1P.log1p(sum) : 1 + sum;
	}


	/* FIXME: The old phyper() code was basically used in ./qhyper.c as well
	 * -----  We need to sync this again!
	*/
	public static double phyper (int x, int NR, int NB, int n,
		       boolean lower_tail, boolean log_p)
	{
	/* Sample of  n balls from  NR red  and	 NB black ones;	 x are red */

	    double d, pd;

	    if (NR < 0 || NB < 0 || Double.isInfinite(NR + NB) || n < 0 || n > NR + NB)
	    	return Double.NaN;

	    if (x * (NR + NB) > n * NR) {
		/* Swap tails.	*/
			int oldNB = NB;
			NB = NR;
			NR = oldNB;
			x = n - x - 1;
			lower_tail = !lower_tail;
	    }

	    if (x < 0) {
	    	return DPQ.R_DT_0(lower_tail, log_p);
	    }
	    
	    if (x >= NR || x >= n) {
	    	return DPQ.R_DT_1(lower_tail, log_p);
	    }

	    d  = dhyper (x, NR, NB, n, log_p);
	    pd = pdhyper(x, NR, NB, n, log_p);

	    return log_p ? DPQ.R_DT_Log(d + pd, lower_tail) : DPQ.R_D_Lval(d * pd, lower_tail);
	}

}