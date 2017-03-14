package org.rfoundation.R.nmath;

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
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
 *    double ftrunc(double x);
 *
 *  DESCRIPTION
 *
 *    Truncation toward zero.
 */

public class Trunc {
	public static double trunc(double x) {
	        if(x >= 0) return Math.floor(x);
	        else return Math.ceil(x);
	}
	
	public static double trunc(int x) {
        if(x >= 0) return Math.floor(x);
        else return Math.ceil(x);
	}
}
