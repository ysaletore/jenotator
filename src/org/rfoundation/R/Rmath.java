package org.rfoundation.R;

/* -*- C -*-
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2011  The R Development Core Team
 *  Copyright (C) 2004       The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *

 * Rmath.h  should contain ALL headers from R's C code in `src/nmath'
   -------  such that ``the Math library'' can be used by simply

   ``#include <Rmath.h> ''

   and nothing else.
*/
public class Rmath {
	/* EPSILONS */
	// DBL_EPSILON = 2^-52 = 2.220446049250313e-16
	public static double DBL_EPSILON = 2.220446049250313e-16;
	
	
	/* ----- The following constants and entry points are part of the R API ---- */

	/* 30 Decimal-place constants */
	/* Computed with bc -l (scale=32; proper round) */

	/* SVID & X/Open Constants */
	/* Names from Solaris math.h */

	public static double M_E = 2.718281828459045235360287471353;	/* e */
	public static double M_LOG2E = 1.442695040888963407359924681002;	/* log2(e) */
	public static double M_LOG10E = 0.434294481903251827651128918917;	/* log10(e) */
	public static double M_LN2 = 0.693147180559945309417232121458;	/* ln(2) */
	public static double M_LN10 = 2.302585092994045684017991454684;	/* ln(10) */
	public static double M_PI = 3.141592653589793238462643383280;	/* pi */
	public static double M_2PI = 6.283185307179586476925286766559;	/* 2*pi */
	public static double M_PI_2 = 1.570796326794896619231321691640;	/* pi/2 */
	public static double M_PI_4 = 0.785398163397448309615660845820;	/* pi/4 */
	public static double M_1_PI = 0.318309886183790671537767526745;	/* 1/pi */
	public static double M_2_PI = 0.636619772367581343075535053490;	/* 2/pi */
	public static double M_2_SQRTPI = 1.128379167095512573896158903122;	/* 2/sqrt(pi) */
	public static double M_SQRT2 = 1.414213562373095048801688724210;	/* sqrt(2) */
	public static double M_SQRT1_2 = 0.707106781186547524400844362105;	/* 1/sqrt(2) */

	/* R-Specific Constants */
	public static double M_SQRT_3 = 1.732050807568877293527446341506;	/* sqrt(3) */
	public static double M_SQRT_32 = 5.656854249492380195206754896838;	/* sqrt(32) */
	public static double M_LOG10_2 = 0.301029995663981195213738894724;	/* log10(2) */
	public static double M_SQRT_PI = 1.772453850905516027298167483341;	/* sqrt(pi) */
	public static double M_1_SQRT_2PI = 0.398942280401432677939946059934;	/* 1/sqrt(2pi) */
	public static double M_SQRT_2dPI = 0.797884560802865355879892119869;	/* sqrt(2/pi) */
	public static double M_LN_SQRT_PI = 0.572364942924700087071713675677;	/* log(sqrt(pi))
		   == log(pi)/2 */

	public static double M_LN_SQRT_2PI = 0.918938533204672741780329736406;	/* log(sqrt(2*pi))
 	 		== log(2*pi)/2 */

	public static double M_LN_SQRT_PId2 = 0.225791352644727432363097614947;	/* log(sqrt(pi/2)) */
}
