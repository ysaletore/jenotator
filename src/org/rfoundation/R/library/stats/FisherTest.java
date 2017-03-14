package org.rfoundation.R.library.stats;

import org.rfoundation.R.nmath.Hyper;

public class FisherTest {
	public static enum Alternative {
		TWO_SIDED,
		LESS,
		GREATER
	}
	
	public static double test(int[][] table) {
		return test(table, Alternative.TWO_SIDED);
	}
	
	public static double test(int[][] table, Alternative alt) {
		if(table.length == 2 && table[0].length == 2 && table[1].length == 2) {
			int m = table[0][0] + table[1][0];
			int n = table[0][1] + table[1][1];
			int k = table[0][0] + table[0][1];
			int x = table[0][0];
			int lo = Math.max(0, k - n);
			int hi = Math.min(k, m);
			double or = 1;
			
			double logdc[] = new double[hi - lo + 1];
			int[] support = new int[hi - lo + 1];
			for(int i = 0; i < logdc.length; i++) {
				support[i] = lo + i;
				logdc[i] = Hyper.dhyper(lo+i, m, n, k, true);
			}
			
			switch(alt) {
				case LESS:
					return pnhyper(x, x, m, n, k, logdc, support, or, false);
				case GREATER:
					return pnhyper(x, x, m, n, k, logdc, support, or, true);
				case TWO_SIDED:
					if(or == 0) {
						return (x == lo) ? 1 : 0;
					} else if(or == Double.POSITIVE_INFINITY) {
						return (x == hi) ? 1 : 0;
					} else {
						double relErr = 1 + Math.pow(10, -7);
						double d[] = dnhyper(logdc, support, or);
						double sum = 0;
						for(int i = 0; i < d.length; i++) {
							if(d[i] <= d[x - lo]) {
								sum += d[i] * relErr;
							}
						}
						
						return Math.max(0, Math.min(sum, 1.0));
					}
				default:
					throw new IllegalArgumentException("Unsupported Alternative test.");
			}
		} else {
			throw new IllegalArgumentException("This test has only been implemented for 2x2 tables.");
		}
	}
	
	private static double[] dnhyper(double[] logdc, int[] support, double ncp) {
		double[] array = new double[support.length];
		
		double max = Double.MIN_VALUE;
		for(int i = 0; i < support.length; i++) {
			array[i] = logdc[i] + Math.log(ncp) * support[i];
			max = Math.max(array[i], max);
		}
		
		for(int i = 0; i < support.length; i++) {
			array[i] = Math.exp(array[i] - max);
		}
		
		double sum = 0.0;
		for(int i = 0; i < support.length; i++) {
			sum += array[i];
		}
		
		for(int i = 0; i < support.length; i++) {
			array[i] /= sum;
		}
		
		return array;
	}
	
	/*
	private static double mnhyper(double[] logdc, int[] support, double ncp) {
		if(ncp == 0) {
			return support[0]; // return lo
		} else if(ncp == Double.POSITIVE_INFINITY) {
			return support[support.length - 1]; // return hi
		} else {
			double sum = 0.0;
			double[] dnhyper_ncp = dnhyper(logdc, support, ncp);
			
			for(int i = 0; i < support.length; i++) {
				sum += support[i] * dnhyper_ncp[i];
			}
			
			return sum;
		}
	}
	*/
	
	private static double pnhyper(double q, int x, int m, int n, int k, 
			double[] logdc, int[] support, double ncp, boolean upper_tail) {
		if(ncp == 1) {
			if(upper_tail) {
				return Hyper.phyper(x - 1, m, n, k, false, false);
			} else {
				return Hyper.phyper(x,     m, n, k, true,  false);
			}
		} else if(ncp == 0) {
			if(upper_tail) {
				return (q <= support[0]) ? 1 : 0;
			} else {
				return (q >= support[0]) ? 1 : 0;
			}
		} else if(ncp == Double.POSITIVE_INFINITY){
			if(upper_tail) {
				return q <= support[support.length - 1] ? 1 : 0;
			} else {
				return q >= support[support.length - 1] ? 1 : 0;
			}
		} else {
			int min = 0, max = support.length - 1;
			if(upper_tail) {
				while(support[min] < q && min < max) {
					min++;
				}
			} else {
				while(support[max] > q && max > min) {
					max--;
				}
			}
			
			double[] dnhyper_ncp = dnhyper(logdc, support, ncp);
			double sum = 0.0;
			if(max >= min) {
				for(int i = min; i <= max; i++) {
					sum += dnhyper_ncp[i];
				}
			}
			
			return sum;
		}
	}
}
