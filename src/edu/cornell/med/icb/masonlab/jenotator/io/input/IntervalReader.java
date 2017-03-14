package edu.cornell.med.icb.masonlab.jenotator.io.input;

import java.util.List;

import edu.cornell.med.icb.masonlab.jenotator.model.interval.Interval;

public interface IntervalReader {
	public boolean hasNext();
	public Interval next();
	public List<Interval> nextSpliced();
}
