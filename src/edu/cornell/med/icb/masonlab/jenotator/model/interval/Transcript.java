package edu.cornell.med.icb.masonlab.jenotator.model.interval;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import edu.cornell.med.icb.masonlab.jenotator.hibernate.ucsc.model.RefGene;

public class Transcript extends Bed6Interval implements StrandedInterval {
	protected List<Interval> introns, cds, utr5, utr3, exons;
	int intronLength, cdsLength, utr5Length, utr3Length;
	protected int cdsStart, cdsEnd, nExons;
	protected String officialName;
	
	public Transcript(RefGene gene) {
		super(gene.getChrom(), gene.getTxStart(), gene.getTxEnd(), 
				gene.getName(), 0, Strand.parseStrand(gene.getStrand()));
		this.cdsStart = gene.getCdsStart();
		this.cdsEnd = gene.getCdsEnd();
		this.officialName = gene.getName2();
		this.nExons = gene.getExonCount();
		
		String[] exonStarts = gene.getExonStarts().split(",");
		String[] exonEnds = gene.getExonEnds().split(",");
		
		this.introns = new ArrayList<Interval>();
		this.cds = new ArrayList<Interval>();
		this.utr5 = new ArrayList<Interval>();
		this.utr3 = new ArrayList<Interval>();
		this.exons = new ArrayList<Interval>();
		
		for(int i = 0; i < this.nExons; i++) {
			Interval exon = new Bed3Interval(chromosome, 
					Integer.parseInt(exonStarts[i]), 
					Integer.parseInt(exonEnds[i]));
			exons.add(exon);
		}
		
		// split up exon into gene body features
		Iterator<Interval> exonIt = exons.iterator();
		Interval exon = exonIt.next();
		while(exon.getEnd() < cdsStart && exonIt.hasNext()) {
			// still in the 5' UTR
			utr5.add(exon);
			utr5Length += exon.getLength();
			exon = exonIt.next();
		}
		
		// now split the exon at the split point
		utr5.add(new Bed3Interval(chromosome, exon.getStart(), cdsStart));
		utr5Length += cdsStart - exon.getStart();
		
		// now add the CDS exons
		while(exon.getEnd() < cdsEnd && exonIt.hasNext()) {
			cds.add(new Bed3Interval(chromosome, Math.max(cdsStart, exon.getStart()), exon.getEnd()));
			cdsLength += exon.getEnd() - Math.max(cdsStart, exon.getStart());
			exon = exonIt.next();
		}
		
		// split the cds exon into cds and 3'UTR
		//TODO: Check cdsEnd OBOB
		cds.add(new Bed3Interval(chromosome, Math.max(cdsStart, exon.getStart()), cdsEnd));
		cdsLength += cdsEnd - Math.max(cdsStart, exon.getStart());
		
		// add the CDS introns into the introns set
		// NOTE: Only CDS introns are considered introns for our purposes
		// introns within the 5' and 3' UTRs are ignored!'
		if(cds.size() > 1) {
			Interval prev = cds.get(0);
			for(int i = 1; i < cds.size(); i++) {
				Interval current = cds.get(i);
				introns.add(new Bed3Interval(chromosome, prev.getEnd(), current.getStart()));
				intronLength += current.getStart() - prev.getEnd();
				prev = current;
			}
		}
		
		utr3.add(new Bed3Interval(chromosome, cdsEnd, exon.getEnd()));
		utr3Length += exon.getEnd() - cdsEnd;
		
		// add the remaining exons into the 3' UTR
		while(exonIt.hasNext()) {
			exon = exonIt.next();
			utr3.add(exon);
			utr3Length += exon.getLength();
		}
		
		/*
		System.out.println("-----------------------");
		System.out.println(this);
		System.out.println("UTR5: " + utr5Length);
		for(Interval interval : utr5) {
			System.out.println("UTR5: " + interval);
		}
		System.out.println("CDS: " + cdsLength);
		for(Interval interval : cds) {
			System.out.println("CDS: " + interval);
		}
		System.out.println("UTR3: " + utr3Length);
		for(Interval interval : utr3) {
			System.out.println("UTR3: " + interval);
		}
		System.out.println("-----------------------");
		*/
	}
	
	public int getIntronLength() {
		return intronLength;
	}

	public int getCDSLength() {
		return cdsLength;
	}

	public int get5UTRLength() {
		return utr5Length;
	}

	public int get3UTRLength() {
		return utr3Length;
	}

	public int getNExons() {
		return this.nExons;
	}

	public List<Interval> getExons() {
		return exons;
	}
	
	public List<Interval> getIntrons() {
		return introns;
	}
	
	public List<Interval> getCDS() {
		return cds;
	}
	
	public List<Interval> get5UTR() {
		return utr5;
	}
	
	public List<Interval> get3UTR() {
		return utr3;
	}

	public int getCdsStart() {
		return cdsStart;
	}

	public void setCdsStart(int cdsStart) {
		this.cdsStart = cdsStart;
	}

	public int getCdsEnd() {
		return cdsEnd;
	}

	public void setCdsEnd(int cdsEnd) {
		this.cdsEnd = cdsEnd;
	}

	public String getOfficialName() {
		return officialName;
	}

	public void setOfficialName(String officialName) {
		this.officialName = officialName;
	}
}
