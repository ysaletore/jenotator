package edu.cornell.med.icb.masonlab.jenotator.annotation.refseq;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.hibernate.Criteria;
import org.hibernate.Session;
import org.hibernate.criterion.Order;

import edu.cornell.med.icb.masonlab.jenotator.hibernate.ucsc.model.RefGene;
import edu.cornell.med.icb.masonlab.jenotator.hibernate.util.UCSCHibernateUtil;

public class GetCurrentRefGeneFromUCSC {
	@SuppressWarnings("unchecked")
	public static List<RefGene> get(String genome) {
		Session session = UCSCHibernateUtil.getSessionFactory(genome).openSession();
		Criteria cr = session.createCriteria(RefGene.class);
		cr.addOrder(Order.asc("chrom"));
		cr.addOrder(Order.asc("txStart"));
		cr.addOrder(Order.asc("name"));
		return (List<RefGene>) cr.list();
	}
	
	public static void main(String[] args) {
		CommandLineParser parser = new GnuParser();
		Options options = new Options();
		buildOptions(options);
		
		try {
			CommandLine cmd = parser.parse(options, args, true);
			String outputFilename = cmd.getOptionValue("output");
			String outputType = cmd.getOptionValue("type");
			String genome = cmd.getOptionValue("genome");
			List<RefGene> refGenes = get(genome);
			
			if(outputType.equalsIgnoreCase("bed")) {
				
			} else if(outputType.equalsIgnoreCase("refGene")) {
				File file = new File(outputFilename);
				PrintStream output = new PrintStream(file);
				for(RefGene refGene : refGenes) {
					output.println(refGene);
				}
				output.close();
			} else if(outputType.equalsIgnoreCase("gtf")) {
				
			} else {
				throw new IllegalArgumentException("Unsupported output type");
			}
		} catch(ParseException e) {
			e.printStackTrace();
			System.err.println("Improper arguments");
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("Jenotator: " + GetCurrentRefGeneFromUCSC.class.getName(), options, true);
			System.exit(1); // so that make stops running
		} catch (IOException e) {
			e.printStackTrace();
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("Jenotator: " + GetCurrentRefGeneFromUCSC.class.getName(), options, true);
			System.exit(1); // so that make stops running
		} catch (Throwable t) {
			t.printStackTrace();
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("Jenotator: " + GetCurrentRefGeneFromUCSC.class.getName(), options, true);
			System.exit(1); // so that make stops running
		}
	}
	
	@SuppressWarnings("static-access")
	public static void buildOptions(Options options) {
		Option genome = OptionBuilder.withArgName("genome")
				  .isRequired(true)
				  .withLongOpt("genome")
				  .withDescription("genome")
				  .hasArg()
				  .create('g');
		options.addOption(genome);
		
		Option output = OptionBuilder.withArgName("output filename")
								  .isRequired(true)
								  .withLongOpt("output")
								  .withDescription("Output filename")
								  .hasArg()
								  .create('o');
		options.addOption(output);
		
		Option type = OptionBuilder.withArgName("output type")
				  .isRequired(true)
				  .withLongOpt("type")
				  .withDescription("Output file type")
				  .hasArg()
				  .create('t');
		options.addOption(type);
	}
}
