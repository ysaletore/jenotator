package edu.cornell.med.icb.masonlab.jenotator.model.interval;

public class Bed3Interval implements Interval {
	protected final String chromosome;
	protected final int start;
	protected final int end;
	
	public Bed3Interval(String chromosome, int start, int end) {
		this.chromosome = chromosome;
		this.start = start;
		this.end = end;
	}
	
	@Override
	public String getChromosome() {
		return chromosome;
	}
	
	@Override
	public int getStart() {
		return start;
	}
	
	@Override
	public int getEnd() {
		return end;
	}
	
	@Override
	public boolean overlaps(Interval o) {
		return this.chromosome.equals(o.getChromosome()) && 
				this.getStart() < o.getEnd() &&
				o.getStart() < this.getEnd();
	}
	
	public boolean isLeft(Interval o) {
		return this.chromosome.equals(o.getChromosome()) &&
				this.getEnd() - 1 < o.getStart();
	}
	
	public boolean isRight(Interval o) {
		return this.chromosome.equals(o.getChromosome()) &&
				this.getStart() >= o.getEnd();
	}
	
	public boolean isLeft(int n) {
		return this.getEnd() - 1 < n;
	}
	
	public boolean isRight(int n) {
		return this.getStart() >= n;
	}
	
	@Override
	public String toString() {
		return this.chromosome + "\t" + this.start + "\t" +
				this.end;
	}
	
	@Override
	public int getLength() {
		return Math.abs(this.getEnd() - this.getStart());
	}
	
	@Override
	public int compareTo(Interval o) {
		if(!this.chromosome.equals(o.getChromosome())) {
			return this.chromosome.compareTo(o.getChromosome());
		}
		
		return ((Integer) this.getStart()).compareTo(o.getStart());
	}
	
	@Override
	public int hashCode() {
		return this.chromosome.hashCode() * 3
				+ ((Integer) this.getStart()).hashCode() * 5
				+ ((Integer) this.getEnd()).hashCode() * 7;
	}
	
	@Override
	public boolean equals(Object other) {
		if(other == null || !(other instanceof Interval)) {
			return false;
		} else {
			Interval o = (Interval) other;
			return this.chromosome.equals(o.getChromosome()) &&
					this.getStart() == o.getStart() &&
					this.getEnd() == o.getEnd();
		}
	}

	@Override
	public int overlapBP(Interval i) {
		if(!this.overlaps(i)) {
			return 0;
		} else {
			int start = Math.max(this.getStart(), i.getStart());
			int end = Math.min(this.getEnd(), i.getEnd());
			return start - end;
		}
	}
}
