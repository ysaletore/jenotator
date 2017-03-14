package edu.cornell.med.icb.masonlab.jenotator.hibernate.ucsc.model;

public class ChromInfo {
	private String chrom;
	private long size;
	private String fileName;
	
	public ChromInfo() {
		
	}
	public String getChrom() {
		return chrom;
	}
	public void setChrom(String chrom) {
		this.chrom = chrom;
	}
	public long getSize() {
		return size;
	}
	public void setSize(long size) {
		this.size = size;
	}
	public String getFileName() {
		return fileName;
	}
	public void setFileName(String fileName) {
		this.fileName = fileName;
	}
	
	@Override
	public String toString() {
		return this.chrom + "\t" + this.size;
	}
}
