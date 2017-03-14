package edu.cornell.med.icb.masonlab.jenotator.collections;

import java.util.Collection;
import java.util.Iterator;
import java.util.NoSuchElementException;

import edu.cornell.med.icb.masonlab.jenotator.model.interval.Interval;

public class AugmentedIntervalTreesByChromosomeIterator implements Iterator<Interval> {
	protected AugmentedIntervalTree tree;
	protected Collection<Iterator<Interval>> iterators;
	protected Iterator<Iterator<Interval>> iteratorIt;
	protected Iterator<Interval> currIt;
	
	public AugmentedIntervalTreesByChromosomeIterator(Collection<Iterator<Interval>> iterators) {
		this.iterators = iterators;
		this.iteratorIt = iterators.iterator();
		this.currIt = null;
	}

	@Override
	public boolean hasNext() {
		while(this.currIt == null || !this.currIt.hasNext()) {
			if(this.iteratorIt.hasNext()) {
				this.currIt = iteratorIt.next();
			} else {
				return false;
			}
		}
		
		return true;
	}

	@Override
	public Interval next() {
		 if (!hasNext()) throw new NoSuchElementException();
		 return currIt.next();
	}

	@Override
	public void remove() {
		this.currIt.remove();
	}
}
