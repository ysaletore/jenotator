package edu.cornell.med.icb.masonlab.jenotator.annotation.profiling.transcriptome;

import java.util.HashMap;
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
public class ProfileCalculator {
	protected boolean plot5UTR;
	protected boolean plotIntrons;
	protected boolean plotCDS;
	protected boolean plot3UTR;
	protected boolean plotFlanking1KB;
	protected CountingMethod countingMethod;
	protected MappingMethod mappingMethod;
	protected TranscriptNormMethod transcriptNormMethod;
	protected int N;
	protected int nBins;
	protected double[][] counters_5UTR;
	protected double[][][] counters_Introns;
	protected double[][][] counters_CDS;
	protected double[][] counters_3UTR;
	protected double[][] counters_flanking1kb5;
	protected double[][] counters_flanking1kb3;
	protected IntervalTree flankingTranscripts5;
	protected IntervalTree flankingTranscripts3;
	protected Map<String, Integer> chrom_sizes;
	protected int[] norm_flanking;
	protected int[] norm_5UTR;
	protected int[][] norm_Introns;
	protected int[][] norm_CDS;
	protected int[] norm_3UTR;
	protected IntervalTree transcripts;
	protected Map<String, Integer> transNorm;
	protected int norm;
	protected int numCounted;
	protected double incrNorm;
		
	public ProfileCalculator(IntervalTree transcripts) {
		this(transcripts, CountingMethod.ALL, MappingMethod.ALL, TranscriptNormMethod.MEAN,
				true, true, false, true, false, 0, null);
	}
	
	public ProfileCalculator(IntervalTree transcripts, boolean plot5UTR, boolean plotCDS, boolean plotIntrons, 
			boolean plot3UTR, boolean plotFlanking1KB, int N, Map<String, Integer> chrom_sizes) {
		this(transcripts, CountingMethod.ALL, MappingMethod.ALL, TranscriptNormMethod.MEAN,
				plot5UTR, plotCDS, plotIntrons, plot3UTR, plotFlanking1KB, N, chrom_sizes);
	}
	
	public ProfileCalculator(IntervalTree transcripts, CountingMethod countingMethod,
			MappingMethod mappingMethod, TranscriptNormMethod transcriptNormMethod,
			boolean plot5UTR, boolean plotCDS, boolean plotIntrons, boolean plot3UTR, boolean plotFlanking1KB, int N, Map<String, Integer> chrom_sizes) {
		this(transcripts, countingMethod, mappingMethod, transcriptNormMethod, plot5UTR, plotCDS, plotIntrons, plot3UTR, plotFlanking1KB, N, 100, chrom_sizes);
	}
	
	public ProfileCalculator(IntervalTree transcripts, CountingMethod countingMethod,
			MappingMethod mappingMethod, TranscriptNormMethod transcriptNormMethod,
			boolean plot5UTR, boolean plotCDS, boolean plotIntrons, boolean plot3UTR, boolean plotFlanking1KB,
			int N, int nBins, Map<String, Integer> chrom_sizes) {
		if(N < 0) {
			throw new IllegalArgumentException("N must be >= 0.");
		}
		this.mappingMethod = mappingMethod;
		this.transcriptNormMethod = transcriptNormMethod;
		this.countingMethod = countingMethod;
		
		this.transcripts = transcripts;
		this.plot5UTR = plot5UTR;
		this.plotIntrons = plotIntrons;
		this.plotCDS = plotCDS;
		this.plot3UTR = plot3UTR;
		this.N = N;
		this.nBins = nBins;
		this.numCounted = 0;
		this.incrNorm = 0;
		this.norm = 0;
		this.plotFlanking1KB = plotFlanking1KB;
		
		this.counters_5UTR = plot5UTR ? new double[N+1][nBins] : null;
		this.counters_Introns = plotIntrons ? new double[N+1][][] : null;
		this.counters_CDS = plotCDS ? new double[N+1][][] : null;
		this.counters_3UTR = plot3UTR ? new double[N+1][nBins] : null;
		
		this.norm_5UTR = plot5UTR ? new int[N+1] : null;
		this.norm_Introns = plotIntrons ? new int[N+1][] : null;
		this.norm_CDS = plotCDS ? new int[N+1][] : null;
		this.norm_3UTR = plot3UTR ? new int[N+1] : null;
		
		if(plotFlanking1KB) {
			this.counters_flanking1kb5 = new double[N+1][nBins];
			this.counters_flanking1kb3 = new double[N+1][nBins];
			this.norm_flanking = new int[N+1];
			this.chrom_sizes = chrom_sizes;
			this.flankingTranscripts5 = new AugmentedIntervalTreesByChromosome();
			this.flankingTranscripts3 = new AugmentedIntervalTreesByChromosome();
		} else {
			this.counters_flanking1kb5 = null;
			this.counters_flanking1kb3 = null;
			this.norm_flanking = null;
			this.flankingTranscripts5 = null;
			this.flankingTranscripts3 = null;
			this.chrom_sizes = null;
		}
		
		if(plotCDS) {
			this.counters_CDS[N] = new double[2*N+1][nBins];
			this.norm_CDS[N] = new int[2*N+1];
			
			for(int i = 0; i < N; i++) {
				this.counters_CDS[i] = new double[2*i+1][nBins];
				this.norm_CDS[i] = new int[2*i+1];
			}
		}
		
		if(plotIntrons) {
			this.counters_Introns[N] = new double[2*N+1][nBins];
			this.norm_Introns[N] = new int[2*N+1];
			
			for(int i = 0; i < N; i++) {
				this.counters_Introns[i] = new double[2*i][nBins];
				this.norm_Introns[i] = new int[2*i];
			}
		}
		
		this.transNorm = new HashMap<String, Integer>();
//		System.out.println("SIZE: " + transcripts.getSize());
		for(Interval interval : transcripts) {
			//System.out.println(interval);
			if(!chrom_sizes.containsKey(interval.getChromosome())) {
				continue; //TODO: Hack to prevent null pointers when ref genome doesn't contain all chromosomes
			}
			String name = ((Transcript) interval).getOfficialName();
			if(!transNorm.containsKey(name)) {
				transNorm.put(name, 1);
			} else {
				transNorm.put(name, transNorm.get(name) + 1);
			}
						
			// use this opportunity to build normalization array
			Transcript transcript = (Transcript) interval;
			int nExons = transcript.getCDS().size();
			int nIntrons = transcript.getIntrons().size();
			int n_exons = (int) Math.ceil((nExons - 1) / 2.0d);
//			int n_introns = (int) Math.ceil((nIntrons - 1) / 2);
			int n = Math.min(n_exons, N);
			
			if(plot5UTR) {
				norm_5UTR[n]++;
				norm++;
			}
			
			if(plot3UTR) {
				norm_3UTR[n]++;
				norm++;
			}
			
			if(plotCDS) {
				for(int i = 0; i < norm_CDS[n].length; i++) {
					norm_CDS[n][i]++;
					norm++;
				}
				
				// fix the overcounting
				if(nExons == norm_CDS[n].length - 1) {
					norm_CDS[n][n]--;
					norm--;
				}
			}
			
			if(plotIntrons) {
				for(int i = 0; i < norm_Introns[n].length; i++) {
					norm_Introns[n][i]++;
					norm++;
				}
				
				// fix the overcounting
				if(nIntrons == norm_Introns[n].length - 1) {
					norm_Introns[n][n]--;
					norm--;
				} else if(nIntrons == norm_Introns[n].length - 2) {
					norm_Introns[n][n]--;
					norm_Introns[n][n+1]--;
					norm--; norm--;
				}
			}
			
			if(plotFlanking1KB) {
				this.flankingTranscripts5.insert(new Bed6Interval(transcript.getChromosome(), 
						Math.max(0, transcript.getStart() - 1000), transcript.getStart(),
						transcript.getOfficialName(), n, transcript.getStrand()
				));
				
				this.flankingTranscripts3.insert(new Bed6Interval(transcript.getChromosome(), 
						transcript.getEnd(), 
						Math.min(this.chrom_sizes.get(transcript.getChromosome()), transcript.getEnd() + 1000),
						transcript.getOfficialName(), n, transcript.getStrand()
				));
				this.norm++;
				this.norm_flanking[n]++;
			}
		}
	}
	
	public void count(List<Interval> intervals) {
		for(Interval interval : intervals) {
			this.count(interval);
		}
	}
	
	public void count(Interval interval) {
		Interval mapped;
		switch(this.mappingMethod) {
			case CENTER:
					mapped = new Bed3Interval(interval.getChromosome(),
							(int) Math.floor((interval.getStart() + interval.getEnd()) / 2),
							(int) Math.floor((interval.getStart() + interval.getEnd()) / 2 + 1));
				break;
			case LEFT:
					mapped = new Bed3Interval(interval.getChromosome(),
							interval.getStart(), interval.getStart() + 1);
				break;
			case RIGHT:
					mapped = new Bed3Interval(interval.getChromosome(),
							interval.getEnd() - 1, interval.getEnd());
				break;
			default:
			case ALL:
				mapped = interval;
				break;
		}
		
		this.numCounted += mapped.getLength() + 1;
		Set<Interval> overlaps = this.transcripts.findOverlaps(mapped);
		if(overlaps.size() > 0) {
			switch(this.countingMethod) {
				case ALL:
					for(Interval trans : overlaps) {
						count(mapped, (Transcript) trans);
					}
					break;
//				case RANDOM:
//					int n = (int) Math.random() * overlaps.size();
//					count(mapped, (Transcript) overlaps.get(n));
//					break;
			}
		}
		
		if(this.plotFlanking1KB) {
			Set<Interval> overlaps5 = this.flankingTranscripts5.findOverlaps(mapped);
			if(overlaps5.size() > 0) {
				switch(this.countingMethod) {
					case ALL:
						for(Interval trans : overlaps5) {
							countFlanking(mapped, (Bed6Interval) trans, this.counters_flanking1kb5, this.counters_flanking1kb3);
						}
						break;
//					case RANDOM:
//						int n = (int) Math.random() * overlaps5.size();
//						countFlanking(mapped, (Bed6Interval) overlaps5.get(n), this.counters_flanking1kb5, this.counters_flanking1kb3);
//						break;
				}
			}
			
			Set<Interval> overlaps3 = this.flankingTranscripts3.findOverlaps(mapped);
			if(overlaps3.size() > 0) {
				switch(this.countingMethod) {
					case ALL:
						for(Interval trans : overlaps3) {
							countFlanking(mapped, (Bed6Interval) trans, this.counters_flanking1kb3, this.counters_flanking1kb5);
						}
						break;
//					case RANDOM:
//						int n = (int) Math.random() * overlaps3.size();
//						countFlanking(mapped, (Bed6Interval) overlaps3.get(n), this.counters_flanking1kb3, this.counters_flanking1kb5);
//						break;
				}
			}
		}
	}
	
	private void count(Interval interval, Transcript transcript) {
		double incr = this.getIncrement(transcript);
		int n = Math.min((int) Math.ceil((transcript.getCDS().size() - 1) / 2.0d), N);
		
		if(plot5UTR) {
			countUTRs(interval, transcript, counters_5UTR, counters_3UTR, transcript.get5UTR(), n, transcript.get5UTRLength(), incr);
		}
		
		if(plotCDS) {
			count_CDS(interval, transcript, counters_CDS, transcript.getCDS(), n, transcript.getCDSLength(), incr);
		}
		
		if(plotIntrons) {
			count_Introns(interval, transcript, counters_Introns, transcript.getIntrons(), n, transcript.getIntronLength(), incr);
		}
		
		if(plot3UTR) {
			countUTRs(interval, transcript, counters_3UTR, counters_5UTR, transcript.get3UTR(), n, transcript.get3UTRLength(), incr);
		}
	}
	
	private void countFlanking(Interval interval, Bed6Interval flanking, double[][] counters, double[][] alt_counters) {
		double incr = this.getIncrement(flanking.getName());
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
			incrNorm += incr;
			if(flanking.isPositiveStrand()) {
				counters[n][i] += incr;
			} else {
				// we are actually in the other flanking going from the opposite end!
				alt_counters[n][nBins-1 - i] += incr;
			}
		}
	}
	
	private void countUTRs(Interval interval, Transcript transcript, double[][] counters, double[][] alt_counters, List<Interval> exons, int n, int length, double incr) {
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
					incrNorm += incr;
					if(transcript.isPositiveStrand()) {
						counters[n][i] += incr;
					} else {
						// we are actually in the other UTR going from the opposite end!
						alt_counters[n][nBins-1 - i] += incr;
					}
				}
			}
			
			pos += exon.getLength();
		}
	}
	
	private void count_CDS(Interval interval, Transcript transcript, double[][][] counters, 
			List<Interval> exons, int n, int length, double incr) {
		// carry out the n separate boxes first
		for(int i = 0; i < n; i++) {
			// 5' end
			Interval exon = exons.get(i);
			if(exon.overlaps(interval)) {
				incrCounters_CDS_Introns(interval, transcript, exon, counters, n, exon.getLength(), incr, i, 0);
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
					incrCounters_CDS_Introns(interval, transcript, exon, counters, n, length, incr, n, pos);
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
				incrCounters_CDS_Introns(interval, transcript, exon, counters, n, exon.getLength(), incr, 2 * n - i, 0);
			}
		}
	}
	
	private void count_Introns(Interval interval, Transcript transcript, double[][][] counters, 
			List<Interval> exons, int n, int length, double incr) {
		if(n == N) {
			// have enough exons to just call the CDS version of the code
			count_CDS(interval, transcript, counters, exons, n, length, incr);
		} else {
			// carry out the n separate boxes first
			// definitely plot all n 5' UTR boxes
			for(int i = 0; i < n; i++) {
				// 5' end
				Interval exon = exons.get(i);
				if(exon.overlaps(interval)) {
					incrCounters_CDS_Introns(interval, transcript, exon, counters, n, exon.getLength(), incr, i, 0);
				}
			}
			
			// carry out the last (n-1) boxes separately
			// this line is here in the event that we're plotting only the last (n-1) introns
			int n2 = Math.min((int) Math.ceil((exons.size() - 1) / 2), N);
			for(int i = 0; i < Math.min(n, n2); i++) {
				// 3' end
				Interval exon = exons.get(exons.size() - i - 1);
				if(exon.overlaps(interval)) {
					incrCounters_CDS_Introns(interval, transcript, exon, counters, n, exon.getLength(), incr, 2 * n - 1 - i, 0);
				}
			}
		}
	}
	
	private void incrCounters_CDS_Introns(Interval interval, Transcript transcript, Interval exon, double[][][] counters, 
			int n, int length, double incr, int box, int pos) {
		// start of overlap relative to exon start
		int start = Math.max(exon.getStart(), interval.getStart()) - exon.getStart();
		
		// end of overlap relative to exon start
		int end = Math.min(exon.getEnd(), interval.getEnd()) - exon.getStart();
		
		// now compute the bins that they overlap
		int sBin = (int) (1.0d * (start + pos) / length * nBins);
		int eBin = (int) Math.ceil(1.0d * (end + pos) / length * nBins);
		
		// increment the counters
		for(int j = Math.max(0, sBin); j <= Math.min(eBin, nBins-1); j++) {
			incrNorm += incr;
			if(transcript.isPositiveStrand()) {
				counters[n][box][j] += incr;
				
			} else {
				counters[n][counters[n].length-1 - box][nBins-1 - j] += incr;
			}
		}
	}
	
	private double getIncrement(Transcript transcript) {
		return getIncrement(transcript.getOfficialName());
	}
	
	private double getIncrement(String name) {
		double incr;
		switch(this.transcriptNormMethod) {
			case MEAN:
				if(name.length() > 0) {
					incr = 1.0d / transNorm.get(name);
				} else {
					incr = 1;
				}
				break;
			default:
			case SUM:
				incr = 1;
				break;
		}
		
		return incr;
	}
	
	public double[][] getCounters_Flanking1KB5() {
		double[][] counters = new double[N+1][nBins];
		for(int i = 0; i <= N; i++) {
			counters[i] = new double[nBins];
			for(int j = 0; j < nBins; j++) {
				counters[i][j] = 1.0d * (1.0d * counters_flanking1kb5[i][j] / incrNorm) / (1.0d * norm_flanking[i] / norm);
				counters[i][j] = 1.0d * (1.0d * counters_flanking1kb5[i][j] / incrNorm) * 100;
				counters[i][j] = 1.0d * (1.0d * counters_flanking1kb5[i][j] / incrNorm) * 100;
			}
		}
		
		return counters;
	}
	
	public double[][] getCounters_Flanking1KB3() {
		double[][] counters = new double[N+1][nBins];
		for(int i = 0; i <= N; i++) {
			counters[i] = new double[nBins];
			for(int j = 0; j < nBins; j++) {
				counters[i][j] = 1.0d * (1.0d * counters_flanking1kb3[i][j] / incrNorm) / (1.0d * norm_flanking[i] / norm);
				counters[i][j] = 1.0d * (1.0d * counters_flanking1kb3[i][j] / incrNorm) * 100;
			}
		}
		
		return counters;
	}

	public double[][] getCounters_5UTR() {
		double[][] counters = new double[N+1][nBins];
		for(int i = 0; i <= N; i++) {
			counters[i] = new double[nBins];
			for(int j = 0; j < nBins; j++) {
				counters[i][j] = 1.0d * (1.0d * counters_5UTR[i][j] / incrNorm) / (1.0d * norm_5UTR[i] / norm);
				counters[i][j] = 1.0d * (1.0d * counters_5UTR[i][j] / incrNorm) * 100;
			}
		}
		
		return counters;
	}

	public double[][][] getCounters_Introns() {
		double[][][] counters = new double[N+1][][];
		for(int i = 0; i <= N; i++) {
			counters[i] = new double[counters_Introns[i].length][];
			for(int j = 0; j < counters_Introns[i].length; j++) {
				counters[i][j] = new double[nBins];
				for(int k = 0; k < nBins; k++) {
					counters[i][j][k] = 1.0d* (1.0d * counters_Introns[i][j][k] / incrNorm) / (1.0d * norm_Introns[i][j] / norm);
					counters[i][j][k] = 1.0d* (1.0d * counters_Introns[i][j][k] / incrNorm) * 100;
				}
			}
		}
		
		return counters;
	}

	public double[][][] getCounters_CDS() {
		double[][][] counters = new double[N+1][][];
		for(int i = 0; i <= N; i++) {
			counters[i] = new double[counters_CDS[i].length][];
			for(int j = 0; j < counters_CDS[i].length; j++) {
				counters[i][j] = new double[nBins];
				for(int k = 0; k < nBins; k++) {
					counters[i][j][k] = 1.0d * (1.0d * counters_CDS[i][j][k] / incrNorm) / (1.0d * norm_CDS[i][j] / norm);
					counters[i][j][k] = 1.0d * (1.0d * counters_CDS[i][j][k] / incrNorm) * 100;
				}
			}
		}
		
		return counters;
	}

	public double[][] getCounters_3UTR() {
		double[][] counters = new double[N+1][nBins];
		for(int i = 0; i <= N; i++) {
			counters[i] = new double[nBins];
			for(int j = 0; j < nBins; j++) {
				counters[i][j] = 1.0d * (1.0d * counters_3UTR[i][j] / incrNorm) / (1.0d * norm_3UTR[i] / norm);
				counters[i][j] = 1.0d * (1.0d * counters_3UTR[i][j] / incrNorm) * 100;
			}
		}
		
		return counters;
	}

	public int[] getNorm_5UTR() {
		return norm_5UTR;
	}

	public int[][] getNorm_Introns() {
		return norm_Introns;
	}

	public int[][] getNorm_CDS() {
		return norm_CDS;
	}

	public int[] getNorm_3UTR() {
		return norm_3UTR;
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
}
