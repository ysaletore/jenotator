package edu.cornell.med.icb.masonlab.jenotator.io.input;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Bed3Interval;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Interval;

public class SamOrBamFileIntervalReader implements IntervalReader {
	SAMFileReader input;
	Iterator<SAMRecord> iterator;
	
	public SamOrBamFileIntervalReader(String filename) throws IOException {
		this(new File(filename));
	}
	
	public SamOrBamFileIntervalReader(File file) throws IOException {
		this.input = new SAMFileReader(file);
		input.setValidationStringency(ValidationStringency.LENIENT);
		this.iterator = this.input.iterator();
	}
	
	/*
	public SamOrBamFileIntervalReader(InputStream inputstream) throws IOException {
		this.bufferedReader = new BufferedReader(new InputStreamReader(inputstream));
		this.readNext();
	}
	*/
	
	@Override
	public boolean hasNext() {
		return this.iterator.hasNext();
	}

	@Override
	public Interval next() {
		SAMRecord record = this.iterator.next();
		// Picard (SAM) uses 1-based, so subtract 1 from the start
		return new Bed3Interval(record.getReferenceName(), record.getAlignmentStart() - 1, record.getAlignmentEnd());
	}
	
	@Override
	public List<Interval> nextSpliced() {
		List<Interval> list = new ArrayList<Interval>(1);
		list.add(this.next());
		return list;
	}
}
