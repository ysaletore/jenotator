package edu.cornell.med.icb.masonlab.jenotator.collections;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.cornell.med.icb.masonlab.jenotator.collections.BinarySearchTreeIterator.TraversalType;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Interval;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.IntervalTree;

public class AugmentedIntervalTreesByChromosome implements IntervalTree {
	protected Map<String, AugmentedIntervalTree> trees;
	protected int size;
	
	public AugmentedIntervalTreesByChromosome() {
		this.trees = new HashMap<String, AugmentedIntervalTree>();
	}

	@Override
	public void insert(Interval interval) {
		if(!this.trees.containsKey(interval.getChromosome())) {
			this.trees.put(interval.getChromosome(), 
					new AugmentedIntervalTree(interval.getChromosome()));
		}
		
		this.trees.get(interval.getChromosome()).insert(interval);
		this.size++;
	}

	@Override
	public Set<Interval> findOverlaps(Interval interval) {
		if(!this.trees.containsKey(interval.getChromosome())) {
			return new HashSet<Interval>();
		}
		
		return this.trees.get(interval.getChromosome()).findOverlaps(interval);
	}
	
	@Override
	public Set<Interval> findOverlaps(List<Interval> intervals) {
		boolean present = false;
		for(Interval i : intervals) {
			if(this.trees.containsKey(i.getChromosome())) {
				present = true;
			}
		}
		
		Set<Interval> set = new HashSet<Interval>(); 
		if(!present) {
			return set;
		}
		
		for(Interval i : intervals) {
			Set<Interval> overlaps = this.trees.get(i.getChromosome()).findOverlaps(intervals);
			for(Interval j : overlaps) {
				set.add(j);
			}
		}
		
		return set;
	}

	@Override
	public int getSize() {
		return this.size;
	}
	
	@Override
	public Iterator<Interval> iterator() {
		return this.iterator(TraversalType.IN_ORDER);
	}
	
	public Iterator<Interval> iterator(TraversalType type) {
		List<Iterator<Interval>> iterators = new ArrayList<Iterator<Interval>>();
		List<String> keys = new ArrayList<String>(this.trees.keySet());
		Collections.sort(keys);
		for(String key : keys) {
			iterators.add(trees.get(key).iterator());
		}
		
		return new AugmentedIntervalTreesByChromosomeIterator(iterators);
	}
}
