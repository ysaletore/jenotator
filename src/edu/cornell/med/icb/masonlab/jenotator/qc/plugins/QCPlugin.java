package edu.cornell.med.icb.masonlab.jenotator.qc.plugins;

import net.sf.samtools.SAMRecord;

public interface QCPlugin {
	public void process(SAMRecord samrecord, String flowcell, int lane);
	public void exportPNG(String filename);
	public void exportPDF(String filename);
	public void exportHTML(String filename);
}
