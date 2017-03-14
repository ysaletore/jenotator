package edu.cornell.med.icb.masonlab.jenotator.io.input;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Bed3Interval;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Interval;

public class SplicedSamOrBamFileIntervalReader implements IntervalReader {
	SAMFileReader input;
	Iterator<SAMRecord> iterator;
	Queue<SAMRecord> queue;
	
	public SplicedSamOrBamFileIntervalReader(String filename) throws IOException {
		this(new File(filename));
	}
	
	public SplicedSamOrBamFileIntervalReader(File file) throws IOException {
		this.input = new SAMFileReader(file);
		input.setValidationStringency(ValidationStringency.LENIENT);
		this.iterator = this.input.iterator();
		this.queue = new ConcurrentLinkedQueue<SAMRecord>();
		this.process();this.process();this.process();this.process();this.process();
	}
	
	/*
	public SamOrBamFileIntervalReader(InputStream inputstream) throws IOException {
		this.bufferedReader = new BufferedReader(new InputStreamReader(inputstream));
		this.readNext();
	}
	*/
	
	private void process() {
		if(iterator.hasNext()) {
			// cycle through until we find a mapped read
			SAMRecord record = iterator.next();
			while(record.getReadUnmappedFlag() && iterator.hasNext()) {
				record = iterator.next();
			}
			queue.add(record);
		}
	}
	
	@Override
	public boolean hasNext() {
		return !this.queue.isEmpty();
	}

	@Override
	public Interval next() {
		this.process();
		SAMRecord record = queue.poll();
		return new Bed3Interval(record.getReferenceName(), record.getAlignmentStart() - 1, record.getAlignmentEnd());
	}
	
	@Override
	public List<Interval> nextSpliced() {
		this.process();
		List<Interval> list = new ArrayList<Interval>();
		SAMRecord record = queue.poll();
		for(AlignmentBlock block : record.getAlignmentBlocks()) {
			int start = (block.getReferenceStart() - 1);
			int end = (block.getReferenceStart() - 1 + block.getLength());
			list.add(new Bed3Interval(record.getReferenceName(), start, end));
		}
		
		return list;
	}
}
