package edu.cornell.med.icb.masonlab.jenotator.statistics;

import java.util.Arrays;

public class LinearRegression {
	public static double[] lm(double[] y, double[] x) {
		if(x.length != y.length){
			throw new IllegalArgumentException("vectors must be of equal length");
		}
		
		double sumx = 0, sumy = 0, xxbar = 0, yybar = 0, xybar = 0;
		int n = x.length;
		for(int i = 0; i < n; i++) {
			sumx += x[i];
			sumy += y[i];
		}
		
		double xbar = sumx / n;
		double ybar = sumy / n;
		
		for(int i = 0; i < n; i++) {
			xxbar += (x[i] - xbar) * (x[i] - xbar);
			yybar += (y[i] - ybar) * (y[i] - ybar);
			xybar += (x[i] - xbar) * (y[i] - ybar);
		}
		
		double beta1 = xybar / xxbar;
		double beta0 = ybar - beta1 * xbar;
		
		double rss = 0.0;
		double ssr = 0.0;
		for(int i = 0; i < n; i++) {
			double fit = beta1 * x[i] + beta0;
			rss += (fit - y[i]) * (fit - y[i]);
			ssr += (fit - ybar) * (fit - ybar);
		}
		
		double R2 = ssr / yybar;
		
		double[] r = {beta0, beta1, R2, rss, ssr};
		
		return r;
	}
	
	public static double[] lm(int[] x, int[] y) {
		return(lm(x, y));
	}
	
	public static void main(String[] args) {
		double[] y = {-0.07885541, 0.77582888, 1.52738944,  0.79387518, -0.40630898, -0.89013517,  1.46783788,  0.67054677, -0.85232780, -0.46911266, -0.10841727,
				-0.05417602,  0.24762643,  0.34439877,  0.16259025,  0.67419868,  0.32892695, -0.89654761,  0.32482147,  0.90578388};
		double[] x = {-0.57593376,  0.98734617,  0.22767377,  0.90639870,  0.04676967, -1.94662863, -0.70970521, -0.87830304,  2.29966543,  0.68684901,  2.59026611,
				-0.72674508, -0.30689601, -0.25590542, -0.51220727, -0.68554848,  1.34142887, -0.71058640,  0.71380272, -1.08001399};
		
		System.out.println(Arrays.toString(lm(y, x)));
	}
}
