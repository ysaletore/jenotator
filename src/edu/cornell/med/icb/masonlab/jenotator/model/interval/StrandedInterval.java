package edu.cornell.med.icb.masonlab.jenotator.model.interval;

public interface StrandedInterval extends Interval {
	public Strand getStrand();
	public boolean isPositiveStrand();
	public boolean isNegativeStrand();
}
