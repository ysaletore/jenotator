package edu.cornell.med.icb.masonlab.jenotator.model.interval;

public interface Interval extends Comparable<Interval> {
	/**
	 * 
	 * @return chromosome
	 */
	public String getChromosome();
	
	/**
	 * 
	 * @return start [bp]
	 */
	public int getStart();
	
	/**
	 * 
	 * @return end [bp]
	 */
	public int getEnd();
	
	/**
	 * 
	 * @return interval length [bp]
	 */
	public int getLength();
	
	/**
	 * Returns if this interval overlaps with interval i (not strand-specific)
	 * @param i	other interval
	 * @return	true if overlaps
	 */
	public boolean overlaps(Interval i);

	/**
	 * Computes the number of bp overlap between this interval and interval i (not strand-specific)
	 * @param i
	 * @return	0 if no overlap; otherwise, number bp in overlap
	 */
	public int overlapBP(Interval i);
}