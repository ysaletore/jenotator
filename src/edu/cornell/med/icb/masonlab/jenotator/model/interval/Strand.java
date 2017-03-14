package edu.cornell.med.icb.masonlab.jenotator.model.interval;

public enum Strand {
	POSITIVE,
	NEGATIVE;
	
	public static Strand parseStrand(String s) {
		return parseStrand(s.charAt(0));
	}
	
	public static Strand parseStrand(Character c) {
		if(c == '-') {
			return NEGATIVE;
		} else {
			return POSITIVE;
		}
	}
	
	public static Strand parseStrand(boolean b) {
		if(b) {
			return POSITIVE;
		} else {
			return NEGATIVE;
		}
	}
	
	@Override
	public String toString() {
		if(this == Strand.POSITIVE) {
			return "+";
		} else {
			return "-";
		}
	}
}