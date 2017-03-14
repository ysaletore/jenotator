package edu.cornell.med.icb.masonlab.jenotator.io.input;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import edu.cornell.med.icb.masonlab.jenotator.model.interval.Bed6Interval;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Strand;

public class Bed6IntervalReader extends Bed3IntervalReader {
	public Bed6IntervalReader(String filename) throws IOException {
		super(filename);
	}
	
	public Bed6IntervalReader(File file) throws IOException {
		super(file);
	}
	
	public Bed6IntervalReader(InputStream inputstream) throws IOException {
		super(inputstream);
	}

	@Override
	public Bed6Interval next() {
		String[] parts = reader.next();
		return new Bed6Interval(parts[0], Integer.parseInt(parts[1]), Integer.parseInt(parts[2]), 
				parts[3], Double.parseDouble(parts[4]), Strand.parseStrand(parts[5]));
	}
}
