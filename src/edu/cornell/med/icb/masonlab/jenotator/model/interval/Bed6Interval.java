package edu.cornell.med.icb.masonlab.jenotator.model.interval;

public class Bed6Interval extends Bed3Interval implements StrandedInterval {
	protected String name;
	protected double score;
	protected Strand strand;
	
	public String getName() {
		return name;
	}

	public double getScore() {
		return score;
	}
	
	public void setScore(double score) {
		this.score = score;
	}

	public Strand getStrand() {
		return strand;
	}
	
	public Bed6Interval(Interval interval) {
		super(interval.getChromosome(), interval.getStart(), interval.getEnd());
	}
	
	public Bed6Interval(String chromosome, int start, int end) {
		super(chromosome, start, end);
	}
	
	public Bed6Interval(String chromosome, int start, int end, String name) {
		super(chromosome, start, end);
		this.name = name;
	}
	
	public Bed6Interval(String chromosome, int start, int end, String name, double score) {
		this(chromosome, start, end, name);
		this.score = score;
	}
	
	public Bed6Interval(String chromosome, int start, int end, String name, double score, Strand strand) {
		this(chromosome, start, end, name, score);
		this.strand = strand;
	}
	
	@Override
	public int hashCode() {
		return this.chromosome.hashCode() * 3
				+ ((Integer) this.start).hashCode() * 5
				+ ((Integer) this.end).hashCode() * 7 
				+ this.name.hashCode() * 11
				+ ((Double) this.score).hashCode() * 13;
	}
	
	@Override
	public boolean equals(Object other) {
		if(other == null) {
			return false;
		} else if(other instanceof Bed3Interval) {
			return super.equals((Bed3Interval) other);
		} else if(other instanceof Bed6Interval){
			Bed6Interval o = (Bed6Interval) other;
			return this.chromosome.equals(o.chromosome) &&
					this.start == o.start &&
					this.end == o.end &&
					this.name.equals(o.name) &&
					this.score == o.score &&
					this.strand == o.strand;
		} else {
			return false;
		}
	}
	
	@Override
	public String toString() {
		return super.toString() + "\t" +
				this.name + "\t" +
				this.score + "\t" + 
				this.strand;
	}

	@Override
	public boolean isPositiveStrand() {
		return this.getStrand() == Strand.POSITIVE;
	}

	@Override
	public boolean isNegativeStrand() {
		return this.getStrand() == Strand.NEGATIVE;
	}
}
