package edu.cornell.med.icb.masonlab.jenotator.qc.plugins;

import net.sf.samtools.SAMRecord;

public class NVC implements QCPlugin {
	private long[] As;
	private long[] Cs;
	private long[] Gs;
	private long[] Ts;
	private long[] Ns;
	
	final byte A = (byte) 'A';
	final byte C = (byte) 'C';
	final byte G = (byte) 'G';
	final byte T = (byte) 'T';
	final byte N = (byte) 'N';
	
	public NVC() {
		As = new long[0];
		Cs = new long[0];
		Gs = new long[0];
		Ts = new long[0];
		Ns = new long[0];
	}

	@Override
	public void process(SAMRecord samrecord, String flowcell, int lane) {
		if(!samrecord.getMateUnmappedFlag()) {
			final byte[] sequence = samrecord.getReadBases();
			if(As.length < sequence.length) {
				long[] newAs = new long[sequence.length];
				long[] newCs = new long[sequence.length];
				long[] newGs = new long[sequence.length];
				long[] newTs = new long[sequence.length];
				long[] newNs = new long[sequence.length];
				
				for(int i = 0; i < As.length; i++) {
					newAs[i] = As[i];
					newCs[i] = Cs[i];
					newGs[i] = Gs[i];
					newTs[i] = Ts[i];
					newNs[i] = Ns[i];
				}
				
				As = newAs;
				Cs = newCs;
				Gs = newGs;
				Ts = newTs;
			}
			
			if(!samrecord.getReadNegativeStrandFlag()) {
				for(int i = 0; i < sequence.length; i++) {
					switch(sequence[i]) {
						case A: 
							As[i]++;
							break;
						case C: 
							Cs[i]++;
							break;
						case G: 
							Gs[i]++;
							break;
						case T: 
							Ts[i]++;
							break;
						default:
						case N: 
							Ns[i]++;
							break;
					}
				}
			} else {
				for(int i = 0; i < sequence.length; i++) {
					switch(sequence[sequence.length - i]) {
						case A: 
							Ts[i]++;
							break;
						case C: 
							Gs[i]++;
							break;
						case G: 
							Cs[i]++;
							break;
						case T: 
							As[i]++;
							break;
						default:
						case N: 
							Ns[i]++;
							break;
					}
				}
			}
		}
	}

	@Override
	public void exportPNG(String filename) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void exportPDF(String filename) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void exportHTML(String filename) {
		// TODO Auto-generated method stub
		
	}

}
