
package edu.cornell.med.icb.masonlab.jenotator.model.interval;

import java.util.List;
import java.util.Set;


public interface IntervalTree extends Iterable<Interval> {
	public abstract void insert(Interval interval);
	public abstract Set<Interval> findOverlaps(Interval interval);
	public abstract Set<Interval> findOverlaps(List<Interval> interval);
	public int getSize();
}
