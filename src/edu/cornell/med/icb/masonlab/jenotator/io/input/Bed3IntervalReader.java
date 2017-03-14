package edu.cornell.med.icb.masonlab.jenotator.io.input;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import edu.cornell.med.icb.masonlab.jenotator.model.interval.Bed3Interval;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Interval;

public class Bed3IntervalReader implements IntervalReader {
	protected TabbedFileReader reader;

	public Bed3IntervalReader(String filename) throws IOException {
		this.reader = new TabbedFileReader(filename);
	}
	
	public Bed3IntervalReader(File file) throws IOException {
		this.reader = new TabbedFileReader(file);
	}
	
	public Bed3IntervalReader(InputStream inputstream) throws IOException {
		this.reader = new TabbedFileReader(inputstream);
	}

	@Override
	public Bed3Interval next() {
		String[] parts = reader.next();
		return new Bed3Interval(parts[0], Integer.parseInt(parts[1]), Integer.parseInt(parts[2]));
	}

	@Override
	public boolean hasNext() {
		return reader.hasNext();
	}
	
	@Override
	public List<Interval> nextSpliced() {
		List<Interval> list = new ArrayList<Interval>(1);
		list.add(this.next());
		return list;
	}
}
