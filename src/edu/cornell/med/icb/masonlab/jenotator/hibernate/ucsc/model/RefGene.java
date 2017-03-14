package edu.cornell.med.icb.masonlab.jenotator.hibernate.ucsc.model;

import edu.cornell.med.icb.masonlab.jenotator.model.interval.Strand;

public class RefGene {
	private int bin;
	private String name;
	private String chrom;
	private Strand strand;
	private int txStart;
	private int txEnd;
	private int cdsStart;
	private int cdsEnd;
	private int exonCount;
	private String exonStarts;
	private String exonEnds;
	private int score;
	private String name2;
	private String cdsStartStat;
	private String cdsEndStat;
	private String exonFrames;
	
	public RefGene() {		
	}
	
	public RefGene(int bin, String name, String chrom, Strand strand,
			int txStart, int txEnd, int cdsStart, int cdsEnd, int exonCount,
			String exonStarts, String exonEnds, int score, String name2,
			String cdsStartStat, String cdsEndStat, String exonFrames) {
		super();
		this.bin = bin;
		this.name = name;
		this.chrom = chrom;
		this.strand = strand;
		this.txStart = txStart;
		this.txEnd = txEnd;
		this.cdsStart = cdsStart;
		this.cdsEnd = cdsEnd;
		this.exonCount = exonCount;
		this.exonStarts = exonStarts;
		this.exonEnds = exonEnds;
		this.score = score;
		this.name2 = name2;
		this.cdsStartStat = cdsStartStat;
		this.cdsEndStat = cdsEndStat;
		this.exonFrames = exonFrames;
	}
	
	public RefGene(String line) {
		String[] parts = line.split("\t");
		this.bin = Integer.parseInt(parts[0]);
		this.name = parts[1];
		this.chrom = parts[2];
		this.strand = Strand.parseStrand(parts[3]);
		this.txStart = Integer.parseInt(parts[4]);
		this.txEnd = Integer.parseInt(parts[5]);
		this.cdsStart = Integer.parseInt(parts[6]);
		this.cdsEnd = Integer.parseInt(parts[7]);
		this.exonCount = Integer.parseInt(parts[8]);
		this.exonStarts = parts[9];
		this.exonEnds = parts[10];
		this.score = Integer.parseInt(parts[11]);
		this.name2 = parts[12];
		this.cdsStartStat = parts[13];
		this.cdsEndStat = parts[14];
		this.exonFrames = parts[15];
	}

	public int getBin() {
		return bin;
	}
	
	public void setBin(int bin) {
		this.bin = bin;
	}
	
	public String getName() {
		return name;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	public String getChrom() {
		return chrom;
	}
	
	public void setChrom(String chrom) {
		this.chrom = chrom;
	}
	
	public char getStrand() {
		return this.strand.toString().charAt(0);
	}
	
	public void setStrand(char strand) {
		this.strand = Strand.parseStrand(strand);
	}
	
	public int getTxStart() {
		return txStart;
	}
	
	public void setTxStart(int txStart) {
		this.txStart = txStart;
	}
	
	public int getTxEnd() {
		return txEnd;
	}
	
	public void setTxEnd(int txEnd) {
		this.txEnd = txEnd;
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
	
	public int getExonCount() {
		return exonCount;
	}
	
	public void setExonCount(int exonCount) {
		this.exonCount = exonCount;
	}
	
	public String getExonStarts() {
		return exonStarts;
	}
	
	public void setExonStarts(String exonStarts) {
		this.exonStarts = exonStarts;
	}
	
	public String getExonEnds() {
		return exonEnds;
	}
	
	public void setExonEnds(String exonEnds) {
		this.exonEnds = exonEnds;
	}
	
	public int getScore() {
		return score;
	}
	
	public void setScore(int score) {
		this.score = score;
	}
	
	public String getName2() {
		return name2;
	}
	
	public void setName2(String name2) {
		this.name2 = name2;
	}
	
	public String getCdsStartStat() {
		return cdsStartStat;
	}
	
	public void setCdsStartStat(String cdsStartStat) {
		this.cdsStartStat = cdsStartStat;
	}
	
	public String getCdsEndStat() {
		return cdsEndStat;
	}
	
	public void setCdsEndStat(String cdsEndStat) {
		this.cdsEndStat = cdsEndStat;
	}
	
	public String getExonFrames() {
		return exonFrames;
	}
	
	public void setExonFrames(String exonFrames) {
		this.exonFrames = exonFrames;
	}
	
	public String toString() {
		StringBuffer buffer = new StringBuffer();
		buffer.append(getBin() + "\t");
		buffer.append(getName() + "\t");
		buffer.append(getChrom() + "\t");
		buffer.append(getStrand() + "\t");
		buffer.append(getTxStart() + "\t");
		buffer.append(getTxEnd() + "\t");
		buffer.append(getCdsStart() + "\t");
		buffer.append(getCdsEnd() + "\t");
		buffer.append(getExonCount() + "\t");
		buffer.append(getExonStarts() + "\t");
		buffer.append(getExonEnds() + "\t");
		buffer.append(getScore() + "\t");
		buffer.append(getName2() + "\t");
		buffer.append(getCdsStartStat() + "\t");
		buffer.append(getCdsEndStat() + "\t");
		buffer.append(getExonFrames());
		return buffer.toString();
	}
}
