package edu.cornell.med.icb.masonlab.jenotator.annotation.profiling.transcriptome;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.cornell.med.icb.masonlab.jenotator.annotation.profiling.CountingMethod;
import edu.cornell.med.icb.masonlab.jenotator.annotation.profiling.MappingMethod;
import edu.cornell.med.icb.masonlab.jenotator.collections.AugmentedIntervalTreesByChromosome;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Bed3Interval;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Bed6Interval;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Interval;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.IntervalTree;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Transcript;

/**
 * 
 * @author Yogesh Saletore
 * TODO: BUGFIXES: 
 * 		- TODO: Overcounting bins when a read spans 2 exons
 */
public class Profiler {
	protected boolean plot5UTR;
	protected boolean plotIntrons;
	protected boolean plotCDS;
	protected boolean plot3UTR;
	protected boolean plotFlanking1KB;
	protected CountingMethod countingMethod;
	protected MappingMethod mappingMethod;
	protected int N;
	protected int nBins;
	protected Map<String, int[][]> counters_5UTR;
	protected Map<String, int[][][]> counters_Introns;
	protected Map<String, int[][][]> counters_CDS;
	protected Map<String, int[][]> counters_3UTR;
	protected Map<String, int[][]> counters_flanking1kb5;
	protected Map<String, int[][]> counters_flanking1kb3;
	protected IntervalTree flankingTranscripts5;
	protected IntervalTree flankingTranscripts3;
	protected Map<String, Integer> chrom_sizes;
	protected IntervalTree transcripts;
	protected int numCounted;
		
	public Profiler(IntervalTree transcripts) {
		this(transcripts, CountingMethod.ALL, MappingMethod.ALL,
				true, true, false, true, false, 0, null);
	}
	
	public Profiler(IntervalTree transcripts, boolean plot5UTR, boolean plotCDS, boolean plotIntrons, 
			boolean plot3UTR, boolean plotFlanking1KB, int N, Map<String, Integer> chrom_sizes) {
		this(transcripts, CountingMethod.ALL, MappingMethod.ALL,
				plot5UTR, plotCDS, plotIntrons, plot3UTR, plotFlanking1KB, N, chrom_sizes);
	}
	
	public Profiler(IntervalTree transcripts, CountingMethod countingMethod,
			MappingMethod mappingMethod,
			boolean plot5UTR, boolean plotCDS, boolean plotIntrons, boolean plot3UTR, boolean plotFlanking1KB, int N, Map<String, Integer> chrom_sizes) {
		this(transcripts, countingMethod, mappingMethod, plot5UTR, plotCDS, plotIntrons, plot3UTR, plotFlanking1KB, N, 100, chrom_sizes);
	}
	
	public Profiler(IntervalTree transcripts, CountingMethod countingMethod,
			MappingMethod mappingMethod,
			boolean plot5UTR, boolean plotCDS, boolean plotIntrons, boolean plot3UTR, boolean plotFlanking1KB,
			int N, int nBins, Map<String, Integer> chrom_sizes) {
		if(N < 0) {
			throw new IllegalArgumentException("N must be >= 0.");
		}
		this.mappingMethod = mappingMethod;
		this.countingMethod = countingMethod;
		
		this.transcripts = transcripts;
		this.plot5UTR = plot5UTR;
		this.plotIntrons = plotIntrons;
		this.plotCDS = plotCDS;
		this.plot3UTR = plot3UTR;
		this.N = N;
		this.nBins = nBins;
		this.numCounted = 0;
		this.plotFlanking1KB = plotFlanking1KB;
		
		this.counters_5UTR = plot5UTR ? new HashMap<String, int[][]>() : null;
		this.counters_Introns = plotIntrons ? new HashMap<String, int[][][]>() : null;
		this.counters_CDS = plotCDS ? new HashMap<String, int[][][]>() : null;
		this.counters_3UTR = plot3UTR ? new HashMap<String, int[][]>() : null;

		if(plotFlanking1KB) {
			this.counters_flanking1kb5 = new HashMap<String, int[][]>();
			this.counters_flanking1kb3 = new HashMap<String, int[][]>();
			this.chrom_sizes = chrom_sizes;
			this.flankingTranscripts5 = new AugmentedIntervalTreesByChromosome();
			this.flankingTranscripts3 = new AugmentedIntervalTreesByChromosome();
		} else {
			this.counters_flanking1kb5 = null;
			this.counters_flanking1kb3 = null;
			this.flankingTranscripts5 = null;
			this.flankingTranscripts3 = null;
			this.chrom_sizes = null;
		}
		
//		System.out.println("SIZE: " + transcripts.getSize());
		for(Interval interval : transcripts) {
			Transcript transcript = (Transcript) interval;
//			out.println(interval);
			if(!chrom_sizes.containsKey(interval.getChromosome())) {
				continue; //TODO: Hack to prevent null pointers when ref genome doesn't contain all chromosomes
			}
			
			String name = ((Transcript) interval).getName();
			if(plotCDS) {
				this.counters_CDS.put(name, new int [N+1][][]);
				
				for(int i = 0; i <= N; i++) {
					this.counters_CDS.get(name)[i] = new int[2*i+1][nBins];
				}
			}
			
			if(plotIntrons) {
				this.counters_Introns.put(name, new int[N+1][][]);
				this.counters_Introns.get(name)[N] = new int[2*N+1][nBins];
				
				for(int i = 0; i < N; i++) {
					this.counters_Introns.get(name)[i] = new int[2*i][nBins];
				}
			}
			
			if(plot5UTR) {
				counters_5UTR.put(name,  new int[N+1][nBins]);
			}
			
			if(plot3UTR) {
				counters_3UTR.put(name, new int[N+1][nBins]);
			}
			
			if(plotFlanking1KB) {
				this.counters_flanking1kb5.put(name, new int[N+1][nBins]);
				this.counters_flanking1kb3.put(name, new int[N+1][nBins]);
				int nExons = transcript.getCDS().size();
				int n_exons = (int) Math.ceil((nExons - 1) / 2.0d);
				int n = Math.min(n_exons, N);
				this.flankingTranscripts5.insert(new Bed6Interval(transcript.getChromosome(), 
						Math.max(0, transcript.getStart() - 1000), transcript.getStart(),
						name, n, transcript.getStrand()
				));
				
				this.flankingTranscripts3.insert(new Bed6Interval(transcript.getChromosome(), 
						transcript.getEnd(), 
						Math.min(this.chrom_sizes.get(transcript.getChromosome()), transcript.getEnd() + 1000),
						name, n, transcript.getStrand()
				));
			}
		}
	}
	
	public void count(Interval interval) {
		List<Interval> list = new ArrayList<Interval>(1);
		list.add(interval);
		count(list);
		
		/*
//		this.numCounted += mapped.getLength() + 1;
		Set<Interval> overlaps = this.transcripts.findOverlaps(list);
		if(overlaps.size() > 0) {
			for(Interval trans : overlaps) {
				count(list, (Transcript) trans);
			}
		}
		
		if(this.plotFlanking1KB) {
			Set<Interval> overlaps5 = this.flankingTranscripts5.findOverlaps(mapped);
			if(overlaps5.size() > 0) {
				for(Interval trans : overlaps5) {
					countFlanking(mapped, (Bed6Interval) trans, this.counters_flanking1kb5, this.counters_flanking1kb3);
				}
			}
			
			Set<Interval> overlaps3 = this.flankingTranscripts3.findOverlaps(mapped);
			if(overlaps3.size() > 0) {
				for(Interval trans : overlaps3) {
					countFlanking(mapped, (Bed6Interval) trans, this.counters_flanking1kb3, this.counters_flanking1kb5);
				}
			}
		}
		*/
	}
	
	public void count(List<Interval> intervals) {
	
//		this.numCounted += mapped.getLength() + 1;
		Set<Interval> set_overlaps = this.transcripts.findOverlaps(intervals);
		boolean matched = false;
		
		if(set_overlaps.size() > 0) {
			for(Interval trans : set_overlaps) {
				matched	= count(intervals, (Transcript) trans) || matched;
			}
		}
		
		if(this.plotFlanking1KB) {
			Set<Interval> overlaps5 = this.flankingTranscripts5.findOverlaps(intervals);
			if(overlaps5.size() > 0) {
				for(Interval trans : overlaps5) {
					countFlanking(intervals, (Bed6Interval) trans, this.counters_flanking1kb5, this.counters_flanking1kb3);
				}
			}
			
			Set<Interval> overlaps3 = this.flankingTranscripts3.findOverlaps(intervals);
			if(overlaps3.size() > 0) {
				for(Interval trans : overlaps3) {
					countFlanking(intervals, (Bed6Interval) trans, this.counters_flanking1kb3, this.counters_flanking1kb5);
				}
			}
		}
	}
	
	private boolean count(List<Interval> intervals, Transcript transcript) {
		int n = Math.min((int) Math.ceil((transcript.getCDS().size() - 1) / 2.0d), N);
		boolean matched = false;

		if(plotCDS) {
			count_CDS(intervals, transcript, counters_CDS, transcript.getCDS(), n, transcript.getCDSLength());
		}

		if(plot3UTR) {
			countUTRs(intervals, transcript, counters_3UTR, counters_5UTR, transcript.get3UTR(), n, transcript.get3UTRLength());
		}
		
		if(plot5UTR) {
			countUTRs(intervals, transcript, counters_5UTR, counters_3UTR, transcript.get5UTR(), n, transcript.get5UTRLength());
		}
		
		
		if(plotIntrons) {
			count_Introns(intervals, transcript, counters_Introns, transcript.getIntrons(), n, transcript.getIntronLength());
		}
		
		return matched;
	}
	
	private static class Counter1 {
		public Map<String, int[][]> counters;
		public String name;
		public int n;
		public int i;
		
		public Counter1(Map<String, int[][]> counters, String name, int n, int i) {
			this.counters = counters;
			this.name = name;
			this.n = n;
			this.i = i;
		}
		
		@Override
		public boolean equals(Object obj) {
			if (this == obj) return true;
		    if (obj == null) return false;
		    if (getClass() != obj.getClass()) return false;
		    Counter1 o = (Counter1) obj;
			return this.name.equals(o.name) && this.n == o.n && this.i == o.i && this.counters == o.counters;
		}
	}
	
	private boolean countFlanking(List<Interval> intervals, Bed6Interval flanking, Map<String, int[][]> counters, Map<String, int[][]> alt_counters) {
		Set<Counter1> set = new HashSet<Counter1>();
		
		for(Interval interval : intervals) {
			// start of overlap relative to exon start
			int start = Math.max(flanking.getStart(), interval.getStart()) - flanking.getStart();
			
			// end of overlap relative to exon start
			int end = Math.min(flanking.getEnd(), interval.getEnd()) - flanking.getStart();;
			
			// now compute the bins that they overlap
			int sBin = (int) (1.0d * start / flanking.getLength() * nBins);
			int eBin = (int) Math.ceil((1.0d * end / flanking.getLength() * nBins));
			
			// increment the counters
			int n = (int) flanking.getScore();
			for(int i = Math.max(0, sBin); i <= Math.min(eBin, nBins-1); i++) {
				if(flanking.isPositiveStrand()) {
					//counters.get(flanking.getName())[n][i]++;
					set.add(new Counter1(counters, flanking.getName(), n, i));
				} else {
					// we are actually in the other flanking going from the opposite end!
					//alt_counters.get(flanking.getName())[n][nBins-1 - i]++;
					set.add(new Counter1(alt_counters, flanking.getName(), n, nBins-1 -i));
				}
			}
		}
		
		// increment the counters
		for(Counter1 c : set) {
			c.counters.get(c.name)[c.n][c.i]++;
		}
		
		return !set.isEmpty();
	}
	
	private boolean countUTRs(List<Interval> intervals, Transcript transcript, Map<String, int[][]> counters, Map<String, int[][]> alt_counters, List<Interval> exons, int n, int length) {
		Set<Counter1> set = new HashSet<Counter1>();
		
		for(Interval interval : intervals) {
			int pos = 0;
			for(Interval exon : exons) {
				if(exon.overlaps(interval)) {
					// start of overlap relative to exon start
					int start = Math.max(exon.getStart(), interval.getStart()) - exon.getStart();
					
					// end of overlap relative to exon start
					int end = Math.min(exon.getEnd(), interval.getEnd()) - exon.getStart();;
					
					// now compute the bins that they overlap
					int sBin = (int) (1.0d * (start + pos) / length * nBins);
					int eBin = (int) Math.ceil((1.0d * (end + pos) / length * nBins));
					
					// increment the counters
					for(int i = Math.max(0, sBin); i <= Math.min(eBin, nBins-1); i++) {
						if(transcript.isPositiveStrand()) {
							//counters.get(transcript.getName())[n][i]++;
							set.add(new Counter1(counters, transcript.getName(), n, i));
						} else {
							// we are actually in the other UTR going from the opposite end!
							//alt_counters.get(transcript.getName())[n][nBins-1 - i]++;
							set.add(new Counter1(alt_counters, transcript.getName(), n, nBins-1 - i));
						}
					}
				}
				
				pos += exon.getLength();
			}
		}
		
		// increment the counters
		for(Counter1 c : set) {
			c.counters.get(c.name)[c.n][c.i]++;
		}
		
		return !set.isEmpty();
	}
	
	private static class Counter2 {
		public Map<String, int[][][]> counters;
		public String name;
		public int n;
		public int i;
		public int j;
		
		public Counter2(Map<String, int[][][]> counters, String name, int n, int i, int j) {
			this.counters = counters;
			this.name = name;
			this.n = n;
			this.i = i;
			this.j = j;
		}
		
		@Override
		public boolean equals(Object obj) {
			if (this == obj) return true;
		    if (obj == null) return false;
		    if (getClass() != obj.getClass()) return false;
		    Counter2 o = (Counter2) obj;
			return this.name.equals(o.name) && this.n == o.n && this.i == o.i && this.counters == o.counters && this.j == o.j;
		}
	}
	
	private boolean count_CDS(List<Interval> intervals, Transcript transcript, Map<String, int[][][]> counters, 
			List<Interval> exons, int n, int length) {
		Set<Counter2> set = new HashSet<Counter2>();
		
		for(Interval interval : intervals) {
			// carry out the n separate boxes first
			for(int i = 0; i < n; i++) {
				// 5' end
				Interval exon = exons.get(i);
				if(exon.overlaps(interval)) {
					incrCounters_CDS_Introns(interval, transcript, exon, counters, n, exon.getLength(), i, 0, set);
				}
				length -= exon.getLength(); // cut out this exon from middle
				
				// adjust the length for the middle section
				exon = exons.get(exons.size() - i - 1);
				length -= exon.getLength(); // cut out this exon from middle
			}
			
			// now compute for the middle box
			if(exons.size() >= 2 * n + 1) {
				int pos = 0;
				for(int i = n; i < exons.size() - n; i++) {
					Interval exon = exons.get(i);
					if(exon.overlaps(interval)) {
						incrCounters_CDS_Introns(interval, transcript, exon, counters, n, length, n, pos, set);
					}
					pos += exon.getLength();
				}
			}
			
			// finish up the 3' end
			// this line is here in the event that we're plotting only the last (n-1) introns
			int n2 = Math.min((int) Math.ceil((exons.size() - 1) / 2.0), N);
			for(int i = 0; i < Math.min(n, n2); i++) {
				// 3' end
				Interval exon = exons.get(exons.size() - i - 1);
				if(exon.overlaps(interval)) {
					incrCounters_CDS_Introns(interval, transcript, exon, counters, n, exon.getLength(), 2 * n - i, 0, set);
				}
			}
		}
		
		// increment the counters
		for(Counter2 c : set) {
			c.counters.get(c.name)[c.n][c.i][c.j]++;
		}
		
		return !set.isEmpty();
	}
	
	private boolean count_Introns(List<Interval> intervals, Transcript transcript, Map<String, int[][][]> counters, 
			List<Interval> exons, int n, int length) {
		Set<Counter2> set = new HashSet<Counter2>();
				
		for(Interval interval : intervals) {
			if(n == N) {
				// have enough exons to just call the CDS version of the code
				// carry out the n separate boxes first
				for(int i = 0; i < n; i++) {
					// 5' end
					Interval exon = exons.get(i);
					if(exon.overlaps(interval)) {
						incrCounters_CDS_Introns(interval, transcript, exon, counters, n, exon.getLength(), i, 0, set);
					}
					length -= exon.getLength(); // cut out this exon from middle
					
					// adjust the length for the middle section
					exon = exons.get(exons.size() - i - 1);
					length -= exon.getLength(); // cut out this exon from middle
				}
				
				// now compute for the middle box
				if(exons.size() >= 2 * n + 1) {
					int pos = 0;
					for(int i = n; i < exons.size() - n; i++) {
						Interval exon = exons.get(i);
						if(exon.overlaps(interval)) {
							incrCounters_CDS_Introns(interval, transcript, exon, counters, n, length, n, pos, set);
						}
						pos += exon.getLength();
					}
				}
				
				// finish up the 3' end
				// this line is here in the event that we're plotting only the last (n-1) introns
				int n2 = Math.min((int) Math.ceil((exons.size() - 1) / 2.0), N);
				for(int i = 0; i < Math.min(n, n2); i++) {
					// 3' end
					Interval exon = exons.get(exons.size() - i - 1);
					if(exon.overlaps(interval)) {
						incrCounters_CDS_Introns(interval, transcript, exon, counters, n, exon.getLength(), 2 * n - i, 0, set);
					}
				}
			} else {
				// carry out the n separate boxes first
				// definitely plot all n 5' UTR boxes
				for(int i = 0; i < n; i++) {
					// 5' end
					Interval exon = exons.get(i);
					if(exon.overlaps(interval)) {
						incrCounters_CDS_Introns(interval, transcript, exon, counters, n, exon.getLength(), i, 0, set);
					}
				}
				
				// carry out the last (n-1) boxes separately
				// this line is here in the event that we're plotting only the last (n-1) introns
				int n2 = Math.min((int) Math.ceil((exons.size() - 1) / 2), N);
				for(int i = 0; i < Math.min(n, n2); i++) {
					// 3' end
					Interval exon = exons.get(exons.size() - i - 1);
					if(exon.overlaps(interval)) {
						incrCounters_CDS_Introns(interval, transcript, exon, counters, n, exon.getLength(), 2 * n - 1 - i, 0, set);
					}
				}
			}
		}
		
		// increment the counters
		for(Counter2 c : set) {
			c.counters.get(c.name)[c.n][c.i][c.j]++;
		}
		
		return !set.isEmpty();
	}
	
	private void incrCounters_CDS_Introns(Interval interval, Transcript transcript, Interval exon, Map<String, int[][][]> counters, 
			int n, int length, int box, int pos, Set<Counter2> set) {
		// start of overlap relative to exon start
		int start = Math.max(exon.getStart(), interval.getStart()) - exon.getStart();
		
		// end of overlap relative to exon start
		int end = Math.min(exon.getEnd(), interval.getEnd()) - exon.getStart();
		
		// now compute the bins that they overlap
		int sBin = (int) (1.0d * (start + pos) / length * nBins);
		int eBin = (int) Math.ceil(1.0d * (end + pos) / length * nBins);
		
		// increment the counters
		for(int j = Math.max(0, sBin); j <= Math.min(eBin, nBins-1); j++) {
			if(transcript.isPositiveStrand()) {
				// counters.get(transcript.getName())[n][box][j]++;
				set.add(new Counter2(counters, transcript.getName(), n, box, j));
			} else {
				// counters.get(transcript.getNam())[n][counters.get(transcript.getName())[n].length-1 - box][nBins-1 - j]++;
				set.add(new Counter2(counters, transcript.getName(), n, counters.get(transcript.getName())[n].length-1 - box, nBins-1 - j));
			}
		}
	}

	public boolean isPlot5UTR() {
		return plot5UTR;
	}

	public boolean isPlotIntrons() {
		return plotIntrons;
	}

	public boolean isPlotCDS() {
		return plotCDS;
	}

	public boolean isPlot3UTR() {
		return plot3UTR;
	}

	public int getN() {
		return N;
	}

	public Map<String, int[][]> getCounters_5UTR() {
		return counters_5UTR;
	}

	public Map<String, int[][][]> getCounters_Introns() {
		return counters_Introns;
	}

	public Map<String, int[][][]> getCounters_CDS() {
		return counters_CDS;
	}

	public Map<String, int[][]> getCounters_3UTR() {
		return counters_3UTR;
	}

	public Map<String, int[][]> getCounters_flanking1kb5() {
		return counters_flanking1kb5;
	}

	public Map<String, int[][]> getCounters_flanking1kb3() {
		return counters_flanking1kb3;
	}
}
