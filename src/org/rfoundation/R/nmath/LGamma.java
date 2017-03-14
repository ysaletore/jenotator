package org.rfoundation.R.nmath;

import org.rfoundation.R.Rmath;

public class LGamma {
	/*
	 *  Mathlib : A C Library of Special Functions
	 *  Copyright (C) 1998 Ross Ihaka
	 *  Copyright (C) 2000-2001 The R Development Core Team
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
	 *    double lgammafn_sign(double x, int *sgn);
	 *    double lgammafn(double x);
	 *
	 *  DESCRIPTION
	 *
	 *    The function lgammafn computes log|gamma(x)|.  The function
	 *    lgammafn_sign in addition assigns the sign of the gamma function
	 *    to the address in the second argument if this is not NULL.
	 *
	 *  NOTES
	 *
	 *    This routine is a translation into C of a Fortran subroutine
	 *    by W. Fullerton of Los Alamos Scientific Laboratory.
	 *
	 *    The accuracy of this routine compares (very) favourably
	 *    with those of the Sun Microsystems portable mathematical
	 *    library.
	 */

	private static double xmax = 2.5327372760800758e+305;
	private static double dxrel = 1.490116119384765696e-8;
	
	public static double lgammafn_sign(double x, int sgn) {
	    double ans, y, sinpiy;

	    if (sgn != 0) sgn = 1;

	    if (x < 0 && Math.floor(-x) % 2 == 0)
		if (sgn != 0) sgn = -1;

	    if (x <= 0 && x == Trunc.trunc(x)) { /* Negative integer argument */
	    	// ML_ERROR(ME_RANGE, "lgamma");
	    	return Double.POSITIVE_INFINITY;/* +Inf, since lgamma(x) = log|gamma(x)| */
	    }

	    y = Math.abs(x);

	    if (y <= 10)
		return Math.log(Math.abs(Gamma.gammafn(x)));
	    /*
	      ELSE  y = |x| > 10 ---------------------- */

	    if (y > xmax) {
	    	// ML_ERROR(ME_RANGE, "lgamma");
	    	return Double.POSITIVE_INFINITY;
	    }

	    if (x > 0) { /* i.e. y = x > 10 */
		    return Rmath.M_LN_SQRT_2PI + (x - 0.5) * Math.log(x) - x + LGammaCor.lgammacor(x);
	    }
	    /* else: x < -10; y = -x */
	    sinpiy = Math.abs(Math.sin(Rmath.M_PI * y));

	    if (sinpiy == 0) { /* Negative integer argument ===
				  Now UNNECESSARY: caught above */
//		MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
	    	return Double.NaN;
	    }

	    ans = Rmath.M_LN_SQRT_PId2 + (x - 0.5) * Math.log(y) 
	    		- x - Math.log(sinpiy) - LGammaCor.lgammacor(y);
	    
	    if(Math.abs((x - Trunc.trunc(x - 0.5)) * ans / x) < dxrel) {

		/* The answer is less than half precision because
		 * the argument is too near a negative integer. */

	    	// FIXME ML_ERROR(ME_PRECISION, "lgamma");
	    }

	    return ans;
	}

	public static double lgammafn(double x)
	{
	    return lgammafn_sign(x, 0);
	}
}
