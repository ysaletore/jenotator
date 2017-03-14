package org.rfoundation.R.nmath;

import org.rfoundation.R.Rmath;

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2000--2007  R Development Core Team
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
 */

/* Utilities for `dpq' handling (density/probability/quantile) */
public class DPQ {
	public static double RD0(boolean give_log) {
		if(give_log) {
			return Double.NEGATIVE_INFINITY;
		} else {
			return 0;
		}
	}
	
	public static double RD1(boolean give_log) {
		if(give_log) {
			return 0;
		} else {
			return 1;
		}
	}
	
	public static double RD_forceint(double x) {
		return Math.floor(x + 0.5);
	}
	
	public static double R_DT_0(boolean lower_tail, boolean give_log) {
		return lower_tail ? RD0(give_log) : RD1(give_log);
	}
	
	public static double R_DT_1(boolean lower_tail, boolean give_log) {
		return lower_tail ? RD1(give_log) : RD0(give_log);
	}
	
	public static double R_D_Lval(double p, boolean lower_tail) {
		return lower_tail ? p : (0.5 - (p) + 0.5);
	}
	
	public static double R_D_Cval(double p, boolean lower_tail) {
		return lower_tail ? (0.5 - (p) + 0.5) : (p);
	}
	
	public static double R_D_val(double x, boolean give_log) {
		return give_log ? Math.log(x) : x;
	}
	
	public static double R_D_qIv(double p, boolean give_log) {
		return give_log ? Math.exp(p) : p;
	}
	
	public static double R_D_exp(double x, boolean give_log) {
		return give_log ? x : Math.exp(x);
	}
	
	public static double R_D_log(double p, boolean give_log) {
		return give_log ? p : Math.log(p);
	}

	public static double R_D_Clog(double p, boolean give_log) {
		return give_log ? Log1P.log1p(-(p)) : (0.5 - (p) + 0.5);
	}

	/* log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x))) : */
	public static double R_Log1_Exp(double x) {
		return ((x) > -Rmath.M_LN2 ? Math.log(-ExpM1.expm1(x)) : Log1P.log1p(-Math.exp(x)));
	}

	/* log(1-exp(x)):  R_D_LExp(x) == (log1p(- R_D_qIv(x))) but even more stable:*/
	public static double R_D_LExp(double x, boolean give_log) {
		return (give_log ? R_Log1_Exp(x) : Log1P.log1p(-x));
	}

	public static double R_DT_val(double x, boolean lower_tail, boolean give_log) {
		return lower_tail ? R_D_val(x, give_log)  : R_D_Clog(x, give_log);
	}
	
	public static double R_DT_Cval(double x, boolean lower_tail, boolean give_log) {
		return lower_tail ? R_D_Clog(x, give_log) : R_D_val(x, give_log);
	}

	/*public static double R_DT_qIv(p)	R_D_Lval(R_D_qIv(p))		 *  p  in qF ! */
	public static double R_DT_qIv(double p, boolean lower_tail, boolean give_log) {
		return (give_log ? (lower_tail ? Math.exp(p) : - ExpM1.expm1(p)) : R_D_Lval(p, lower_tail));
	}

	/*public static double R_DT_CIv(p)	R_D_Cval(R_D_qIv(p))		 *  1 - p in qF */
	public static double R_DT_CIv(double p, boolean lower_tail, boolean give_log) {
		return (give_log ? (lower_tail ? -ExpM1.expm1(p) : Math.exp(p)) : R_D_Cval(p, lower_tail));
	}

	public static double R_DT_exp(double x, boolean lower_tail, boolean give_log) {
		return R_D_exp(R_D_Lval(x, lower_tail), give_log);		/* exp(x) */
	}
	
	public static double R_DT_Cexp(double x, boolean give_log) {
		return R_D_exp(R_D_Cval(x,give_log),give_log);		/* exp(1 - x) */
	}

	public static double R_DT_log(double p, boolean lower_tail, boolean give_log) {
		return lower_tail ? R_D_log(p, give_log) : R_D_LExp(p, give_log); /* log(p) in qF */
	}
	
	public static double R_DT_Clog(double p, boolean lower_tail, boolean give_log) {
		return lower_tail ? R_D_LExp(p, give_log): R_D_log(p, give_log); /* log(1-p) in qF*/
	}
	
	public static double R_DT_Log(double p, boolean lower_tail) {
		return lower_tail ? (p) : R_Log1_Exp(p);
	}
	
	public static double R_D_fexp(double f, double x, boolean give_log) {
		return give_log ? -0.5*Math.log(f)+(x) : Math.exp(x)/Math.sqrt(f);
	}

}
