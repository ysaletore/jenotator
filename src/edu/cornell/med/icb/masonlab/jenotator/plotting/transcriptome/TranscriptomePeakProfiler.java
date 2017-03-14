package edu.cornell.med.icb.masonlab.jenotator.plotting.transcriptome;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import edu.cornell.med.icb.masonlab.jenotator.activity.GetChromosomeSizes;
import edu.cornell.med.icb.masonlab.jenotator.annotation.profiling.CountingMethod;
import edu.cornell.med.icb.masonlab.jenotator.annotation.profiling.MappingMethod;
import edu.cornell.med.icb.masonlab.jenotator.annotation.profiling.transcriptome.Profiler;
import edu.cornell.med.icb.masonlab.jenotator.annotation.profiling.transcriptome.TranscriptNormMethod;
import edu.cornell.med.icb.masonlab.jenotator.annotation.refseq.GetRefGeneFromFile;
import edu.cornell.med.icb.masonlab.jenotator.collections.AugmentedIntervalTreesByChromosome;
import edu.cornell.med.icb.masonlab.jenotator.hibernate.ucsc.model.RefGene;
import edu.cornell.med.icb.masonlab.jenotator.io.input.Bed3IntervalReader;
import edu.cornell.med.icb.masonlab.jenotator.io.input.IntervalReader;
import edu.cornell.med.icb.masonlab.jenotator.io.input.SamOrBamFileIntervalReader;
import edu.cornell.med.icb.masonlab.jenotator.io.input.SplicedSamOrBamFileIntervalReader;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Interval;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.IntervalTree;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Transcript;

public class TranscriptomePeakProfiler {
	private static class SimplePrinter implements Callable<Integer> {
		protected PrintStream ps;
		protected Map<String, int[][]> c;
		protected List<String> names;
		protected int N;
		protected int num_bins;
		
		public SimplePrinter(PrintStream ps2, Map<String, int[][]> c2, List<String> names2, int num_bins2, int N2) {
			this.ps = ps2;
			this.c = c2;
			this.names = names2;
			this.num_bins = num_bins2;
			this.N = N2;
		}
		
		public Integer call() {
			// print header
			for(int j = 0; j <= N; j++) {
				for(String name : names) {
					ps.print(name + "\t");
				}
			}
			ps.println();
			
			for(int bin = 0; bin < num_bins; bin++) {
				for(int j = 0; j <= N; j++) {
					for(String name : names) {
						ps.print(c.get(name)[j][bin] + "\t");
					}
				}
				ps.println();
			}
			
			return new Integer(1);
		}
	}
	
	private static class CDSPrinter implements Callable<Integer> {
		protected PrintStream ps;
		protected Map<String, int[][][]> c;
		protected List<String> names;
		protected int N;
		protected int num_bins;
		
		public CDSPrinter(PrintStream ps2, Map<String, int[][][]> c2, List<String> names2, int num_bins2, int N2) {
			this.ps = ps2;
			this.c = c2;
			this.names = names2;
			this.num_bins = num_bins2;
			this.N = N2;
		}
		
		public Integer call() {
			// print header
			for(int j = 0; j <= N; j++) {
				for(int k = 0; k < 2*j + 1; k++) {
					for(String name : names) {
						ps.print(name + "\t");
					}
				}
			}
			ps.println();
			
			
			for(int bin = 0; bin < num_bins; bin++) {
				for(int j = 0; j <= N; j++) {
					for(int k = 0; k < 2 * j + 1; k++) {
						for(String name : names) {
							ps.print(c.get(name)[j][k][bin] + "\t");
						}
					}
				}
				ps.println();
			}
			
			return new Integer(1);
		}
	}
	
	private static class IntronsPrinter implements Callable<Integer> {
		protected PrintStream ps;
		protected Map<String, int[][][]> c;
		protected List<String> names;
		protected int N;
		protected int num_bins;
		
		public IntronsPrinter(PrintStream ps2, Map<String, int[][][]> c2, List<String> names2, int num_bins2, int N2) {
			this.ps = ps2;
			this.c = c2;
			this.names = names2;
			this.num_bins = num_bins2;
			this.N = N2;
		}
		
		public Integer call() {
			// print header
			for(int j = 0; j <= N; j++) {
				for(int k = 0; k < 2 * j + 1; k++) {
					for(String name : names) {
						ps.print(name + "\t");
					}
				}
			}
			ps.println();
			
			
			for(int bin = 0; bin < num_bins; bin++) {
				for(int j = 0; j < N; j++) {
					for(int k = 0; k < 2 * j; k++) {
						for(String name : names) {
							ps.print(c.get(name)[j][k][bin] + "\t");
						}
					}
				}
				
				for(int k = 0; k < 2 * N + 1; k++) {
					for(String name : names) {
						ps.print(c.get(name)[N][k][bin] + "\t");
					}
				}
				ps.println();
			}
			
			return new Integer(1);
		}
	}
	
	private static class TranscriptomeProfilePlotThread implements Callable<Integer> {
		protected PrintStream ps_utr5;
		protected PrintStream ps_utr3;
		protected PrintStream ps_introns;
		protected PrintStream ps_cds;
		protected PrintStream ps_flanking5;
		protected PrintStream ps_flanking3;
		protected IntervalReader reader;
		protected Profiler profiler;
		protected IntervalTree transcripts;
		protected int num_bins;
		protected int N;
		protected Map<String, Integer> chrom_sizes;
		protected ExecutorService threadPool;
		
		public TranscriptomeProfilePlotThread(boolean plot5UTR, boolean plotCDS, boolean plotIntrons, boolean plot3UTR,boolean plot1KBUp, boolean plot1KBDown, int N,
				int num_bins, String inputFilename, String inputType, String prefix, String transcriptsFilename, String transcriptsFileType, Map<String, Integer> chrom_sizes, ExecutorService threadPool) throws IOException {
			this(plot5UTR, plotCDS, plotIntrons, plot3UTR, plot1KBUp, plot1KBDown, N, num_bins, inputFilename, inputType, prefix, chrom_sizes, threadPool);
			this.transcripts = getTranscripts(transcriptsFilename, transcriptsFileType);
		}
		
		public TranscriptomeProfilePlotThread(boolean plot5UTR, boolean plotCDS, boolean plotIntrons, boolean plot3UTR, boolean plot1KBUp, boolean plot1KBDown, int N,
				int num_bins, String inputFilename, String inputType, String prefix, IntervalTree transcripts, Map<String, Integer> chrom_sizes, ExecutorService threadPool) throws IOException {
			this(plot5UTR, plotCDS, plotIntrons, plot3UTR, plot1KBUp, plot1KBDown, N, num_bins, inputFilename, inputType, prefix, chrom_sizes, threadPool);
			this.transcripts = transcripts;
		}
		
		private TranscriptomeProfilePlotThread(boolean plot5UTR, boolean plotCDS, boolean plotIntrons, boolean plot3UTR, boolean plot1KBUp, boolean plot1KBDown, int N,
				int num_bins, String inputFilename, String inputType, String prefix, Map<String, Integer> chrom_sizes, ExecutorService threadPool) throws IOException {
			if(inputType.equalsIgnoreCase("bed")) {
				reader = new Bed3IntervalReader(inputFilename);
			} else if(inputType.equalsIgnoreCase("sam")) {
				reader = new SplicedSamOrBamFileIntervalReader(inputFilename);
			} else {
				throw new IllegalArgumentException(inputType + " is not a supported input format.");
			}
			
			ps_utr5 = plot5UTR ? new PrintStream(new File(prefix + "5utr.txt")) : null;
			ps_utr3 = plot3UTR ? new PrintStream(new File(prefix + "3utr.txt")) : null;
			ps_introns = plotIntrons ? new PrintStream(new File(prefix + "introns.txt")) : null;
			ps_cds = plotCDS ? new PrintStream(new File(prefix + "cds.txt")) : null;
			ps_flanking5 = plot1KBUp ? new PrintStream(new File(prefix + "flanking5.txt")) : null;
			ps_flanking3 = plot1KBDown ? new PrintStream(new File(prefix + "flanking3.txt")) : null;
			
			this.num_bins = num_bins;
			this.N = N;
			this.chrom_sizes = chrom_sizes;
			this.threadPool = threadPool;
		}
		
		@Override
		public Integer call() {
			try {
				System.out.println("Started");
				// get sorted list of transcript names
				List<String> names = new ArrayList<String>();
				for(Interval transcript : transcripts) {
					if(!chrom_sizes.containsKey(transcript.getChromosome())) {
						continue; //TODO: Hack to prevent null pointers when ref genome doesn't contain all chromosomes
					}
					names.add(((Transcript) transcript).getName());
				}
				Collections.sort(names);
				List<Future<Integer>> futures = new ArrayList<Future<Integer>>();
				try {
					profiler = new Profiler(transcripts, CountingMethod.ALL, MappingMethod.ALL, 
							ps_utr5 != null, ps_cds != null, ps_introns != null, ps_utr3 != null, 
							ps_flanking5 != null || ps_flanking3 != null, N, num_bins, chrom_sizes);
					System.out.println("Running");
					List<Interval> read;
					int counter = 0;
				 	while(reader.hasNext()) {
						read = reader.nextSpliced();
						counter++;
						if(read != null && !read.isEmpty()) {
							profiler.count(read);
						}
						
						if(counter % 1e7 == 0) {
							System.out.println("Processed " + counter + " reads");
						}
					}
				} catch(Exception e) {
					e.printStackTrace();
					throw e;
				}
				
				System.out.println("printing");
				if(ps_flanking5 != null) {
					futures.add(threadPool.submit(new SimplePrinter(ps_flanking5, profiler.getCounters_flanking1kb5(), names, num_bins, N)));
				}
				
				if(ps_utr5 != null) {
					futures.add(threadPool.submit(new SimplePrinter(ps_utr5, profiler.getCounters_5UTR(), names, num_bins, N)));
				}
				
				if(ps_cds != null) {
					futures.add(threadPool.submit(new CDSPrinter(ps_cds, profiler.getCounters_CDS(), names, num_bins, N)));
				}
				
				if(ps_introns != null) {
					futures.add(threadPool.submit(new IntronsPrinter(ps_introns, profiler.getCounters_Introns(), names, num_bins, N)));
				}
				
				if(ps_flanking3 != null) {
					futures.add(threadPool.submit(new SimplePrinter(ps_flanking3, profiler.getCounters_flanking1kb3(), names, num_bins, N)));
				}
				
				if(ps_utr3 != null) {
					futures.add(threadPool.submit(new SimplePrinter(ps_utr3, profiler.getCounters_3UTR(), names, num_bins, N)));
				}
				
				Integer sum = 0;
				for(Future<Integer> future : futures) {
					try {
						sum += future.get();
					} catch(Exception e) {
						e.printStackTrace();
					}
				}
				
				return sum;
			} catch(Throwable T) {
				T.printStackTrace();
			}
			
			return new Integer(-1);
		}
	}
	
	public static void main(String[] args) {
		CommandLineParser parser = new GnuParser();
		Options options = new Options();
		buildOptions(options);
		
		try {
			CommandLine cmd = parser.parse(options, args, true);
			
			int n_samples;
			int n_threads = 1;
			int num_bins = 100;
			int rolling_mean = 5;
			
			boolean plot1KBUp = cmd.hasOption("plot1KBUp");
			boolean plot1KBDown = cmd.hasOption("plot1KBDown");
			boolean plot5UTR = cmd.hasOption("plot5UTR");
			boolean plotCDS = cmd.hasOption("plotCDS");
			boolean plotIntrons = cmd.hasOption("plotIntrons");
			boolean plot3UTR = cmd.hasOption("plot3UTR");
			
			String[] inputFilenames = cmd.getOptionValues("input");
			String[] inputTypes = cmd.getOptionValues("input-type");
			String[] prefixes = cmd.getOptionValues("prefix");
			
			n_samples = inputFilenames.length;
			if(n_samples != inputTypes.length) {
				throw new IllegalArgumentException("# of inputted file types must equal # of inputted files.");
			}
			if(n_samples != prefixes.length) {
				throw new IllegalArgumentException("# of inputted prefixes must equal # of inputted files.");
			}
			
			String[] transcriptsFilenames = cmd.getOptionValues("transcripts");
			String[] transcriptsFileTypes = cmd.getOptionValues("transcripts-type");
			IntervalTree transcripts = null;
			
			if(transcriptsFilenames.length != transcriptsFileTypes.length) {
				throw new IllegalArgumentException("# of transcript filenames must equal # of transcript types.");
			}
			
			String[] genomeSizesFilenames = cmd.getOptionValues("genome-sizes");
			List<Map<String, Integer>> chrom_sizes_list = new ArrayList<Map<String, Integer>>(n_samples);
			if(genomeSizesFilenames.length == 1) {
				Map<String, Integer> chrom_sizes = GetChromosomeSizes.get(genomeSizesFilenames[0]);
				for(int i = 0; i < n_samples; i++) {
					chrom_sizes_list.add(chrom_sizes);
				}
			} else if(genomeSizesFilenames.length == n_samples) {
				for(int i = 0; i < n_samples; i++) {
					chrom_sizes_list.add(GetChromosomeSizes.get(genomeSizesFilenames[i]));
				}
			} else {
				throw new IllegalArgumentException("Inputted # of genome sizes files must be 1 or equal to # of inputted files.");
			}
			
			if(cmd.hasOption("num-threads")) {
				n_threads = Integer.parseInt(cmd.getOptionValue("num-threads"));
			}
			
			if(cmd.hasOption("num-bins")) {
				num_bins = Integer.parseInt(cmd.getOptionValue("num-bins"));
				if(num_bins < 1) {
					throw new IllegalArgumentException("# of bins must be >= 1.");
				}
			}
			
			if(cmd.hasOption("rolling-mean")) {
				rolling_mean = Integer.parseInt(cmd.getOptionValue("rolling-mean"));
				if(rolling_mean < 1) {
					throw new IllegalArgumentException("rolling_mean must be >= 1.");
				}
			}
			
			int N = 0;
			if(cmd.hasOption("N")) {
				N = Integer.parseInt(cmd.getOptionValue("N"));
				if(N < 0) {
					throw new IllegalArgumentException("N must be >= 0.");
				}
			}
			
			String R_col = cmd.getOptionValue("R-col");
			if(R_col == null) {
				R_col = "blue";
			}
			
			String R_lty = cmd.getOptionValue("R-lty");
			if(R_lty == null) {
				R_lty = "solid,dashed";
			}
			
			String R_lwd = cmd.getOptionValue("R-lwd");
			if(R_lwd == null) {
				R_lwd = "1.5";
			}
			
			if(transcriptsFilenames.length == 1) {
				transcripts = getTranscripts(transcriptsFilenames[0], transcriptsFileTypes[0]);
			}
			
			List<String> cdsFilesList = new ArrayList<String>();
			List<String> utr5FilesList = new ArrayList<String>();
			List<String> utr3FilesList = new ArrayList<String>();
			List<String> intronsFilesList = new ArrayList<String>();
			List<String> flanking5FilesList = new ArrayList<String>();
			List<String> flanking3FilesList = new ArrayList<String>();
			
			ExecutorService threadPool = Executors.newFixedThreadPool(n_threads);
			List<Future<Integer>> futures = new ArrayList<Future<Integer>>();
			for(int i = 0; i < inputFilenames.length; i++) {
				String inputType = inputTypes[i];
				String inputFilename = inputFilenames[i];
				
				if(transcripts == null) {
					futures.add(threadPool.submit(new TranscriptomeProfilePlotThread(plot5UTR, plotCDS, plotIntrons, plot3UTR, plot1KBUp, plot1KBDown, N, 
							num_bins, inputFilename, inputType, prefixes[i], transcriptsFilenames[i], transcriptsFileTypes[i], chrom_sizes_list.get(i), threadPool)));
				} else {
					futures.add(threadPool.submit(new TranscriptomeProfilePlotThread(plot5UTR, plotCDS, plotIntrons, plot3UTR, plot1KBUp, plot1KBDown, N, 
							num_bins, inputFilename, inputType, prefixes[i], transcripts, chrom_sizes_list.get(i), threadPool)));
				}
				
				if(plotCDS) cdsFilesList.add(prefixes[i] + "cds.txt");
				if(plot5UTR) utr5FilesList.add(prefixes[i] + "5utr.txt");
				if(plot3UTR) utr3FilesList.add(prefixes[i] + "3utr.txt");
				if(plotIntrons) intronsFilesList.add(prefixes[i] + "introns.txt");
				if(plot1KBUp) flanking5FilesList.add(prefixes[i] + "flanking5.txt");
				if(plot1KBDown) flanking3FilesList.add(prefixes[i] + "flanking3.txt");	
			}
			try {
				Thread.sleep(5*1000);
			} catch(InterruptedException e) {
				e.printStackTrace();
			}
			
			System.out.println("Waiting for threads to complete...");
			for(Future<Integer> future : futures) {
				try {
					Integer i = future.get();
					System.out.println("Thread completed: " + i);
				} catch(Exception e) {
					e.printStackTrace();
				}
			}
			
			threadPool.shutdown();
			
			while(!threadPool.isShutdown() || !threadPool.isTerminated()) {
				System.out.println("Waiting for threads to complete...");
				try {
					Thread.sleep(5000);
				} catch(InterruptedException e) {
					e.printStackTrace();
				}
			}
		} catch(ParseException e) {
			e.printStackTrace();
			System.err.println("Improper arguments");
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("Jenotator: " + TranscriptomeProfilePlot.class.getName(), options, true);
			System.exit(1); // so that make stops running
		} catch(IllegalArgumentException e) {
			e.printStackTrace();
			System.err.println("Improper arguments");
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("Jenotator: " + TranscriptomeProfilePlot.class.getName(), options, true);
			System.exit(1); // so that make stops running
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1); // so that make stops running
		} catch (Throwable t) {
			t.printStackTrace();
			System.exit(1); // so that make stops running
		}
	}
	
	private static IntervalTree getTranscripts(String filename, String type) throws IOException {
		List<RefGene> refGenes = GetRefGeneFromFile.load(filename, type);
		IntervalTree tree = new AugmentedIntervalTreesByChromosome();
		for(RefGene refGene : refGenes) {
			if(refGene.getName().startsWith("NM")) {
				tree.insert(new Transcript(refGene));
			}
		}
		
		return tree;
	}
	
	@SuppressWarnings("static-access")
	public static void buildOptions(Options options) {
		Option input = OptionBuilder.withArgName("input filename")
							  .isRequired(true)
							  .withLongOpt("input")
							  .withDescription("Input filename")
							  .hasArgs()
							  .withValueSeparator(',')
							  .create('i');
		options.addOption(input);
		
		Option N = OptionBuilder.withArgName("N")
				  .withLongOpt("N")
				  .withDescription("# of exons to separate")
				  .hasArg()
				  .create('N');
		options.addOption(N);
		
		Option plotFlankingUp = OptionBuilder.isRequired(false)
				  .withLongOpt("plot1KBUp")
				  .withDescription("Plot 1 KB upstream")
				  .hasArg(false)
				  .create();
		options.addOption(plotFlankingUp);
		
		Option plotFlankingUpDown = OptionBuilder.isRequired(false)
				  .withLongOpt("plot1KBDown")
				  .withDescription("Plot 1 KB downstream")
				  .hasArg(false)
				  .create();
		options.addOption(plotFlankingUpDown);
		
		Option plot5UTR = OptionBuilder.isRequired(false)
				  .withLongOpt("plot5UTR")
				  .withDescription("Plot 5' UTR")
				  .hasArg(false)
				  .create();
		options.addOption(plot5UTR);
		
		Option plotCDS = OptionBuilder.isRequired(false)
				  .withLongOpt("plotCDS")
				  .withDescription("Plot CDS")
				  .hasArg(false)
				  .create();
		options.addOption(plotCDS);
		
		Option plotIntrons = OptionBuilder.isRequired(false)
				  .withLongOpt("plotIntrons")
				  .withDescription("Plot Introns")
				  .hasArg(false)
				  .create();
		options.addOption(plotIntrons);
		
		Option plot3UTR = OptionBuilder.isRequired(false)
				  .withLongOpt("plot3UTR")
				  .withDescription("Plot 3' UTR")
				  .hasArg(false)
				  .create();
		options.addOption(plot3UTR);
		
		Option genome = OptionBuilder.withArgName("genome.sizes file(s)")
				.hasArgs()
				.withValueSeparator(',')
				.isRequired()
				.withLongOpt("genome-sizes")
				.withDescription("Genome Chromosome Sizes - must be 1 argument or length equal to # of samples")
				.create();
		options.addOption(genome);
		
		Option type = OptionBuilder.withArgName("input type(s)")
				  .isRequired(true)
				  .withLongOpt("input-type")
				  .withDescription("input file type - must be of length equal to # of samples")
				  .hasArgs()
				  .withValueSeparator(',')
				  .create();
		options.addOption(type);
		
		Option prefix = OptionBuilder.withArgName("prefixes")
				.hasArgs()
				.withValueSeparator(',')
				.isRequired(true)
				.withLongOpt("prefix")
				.withDescription("file prefix(s) - must be of length equal to # of samples")
				.create();
		options.addOption(prefix);
		
		Option transcripts = OptionBuilder.withArgName("transcripts file(s)")
				  .isRequired(true)
				  .withLongOpt("transcripts")
				  .withDescription("transcripts - must be of length 1 or equal to # of samples")
				  .hasArgs()
				  .withValueSeparator(',')
				  .create();
		options.addOption(transcripts);
		
		Option type2 = OptionBuilder.withArgName("transcripts type(s)")
				  .isRequired(true)
				  .withLongOpt("transcripts-type")
				  .withDescription("transcripts type - must be list of length equal to number of inputted transcript files.")
				  .hasArgs()
				  .withValueSeparator(',')
				  .create();
		options.addOption(type2);

		Option nthreads = OptionBuilder.withArgName("# threads")
				  .isRequired(false)
				  .withLongOpt("num-threads")
				  .withDescription("# of threads to use - useful if processing more than one sample.")
				  .hasArg()
				  .create('t');
		options.addOption(nthreads);
		
		Option nbins = OptionBuilder.withArgName("# bins")
				  .isRequired(false)
				  .withLongOpt("num-bins")
				  .withDescription("# of bins to use")
				  .hasArg()
				  .create();
		options.addOption(nbins);
	}
}
