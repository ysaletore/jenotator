package edu.cornell.med.icb.masonlab.jenotator.collections;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.cornell.med.icb.masonlab.jenotator.model.interval.Interval;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.IntervalTree;

public class AugmentedIntervalTree extends AVLTree<Interval> implements IntervalTree  {
	protected String chromosome;
	
	public AugmentedIntervalTree(String chromosome) {
		this.chromosome = chromosome;
	}
	
	@Override
	public void insert(Interval interval) {
		if(!this.chromosome.equals(interval.getChromosome())) {
			throw new ChromosomeMismatchException();
		}
		
		this.root = insert((AugmentedIntervalTreeNode) this.root, interval);
		this.size++;
	}
	
	/** 
	 * insert: 					insert interval s into tree
	 * @param 	AVLNode	node	node to start insertion at
	 * @param	String	s		string to insert
	 * @return	AVLNode			the new node
	 */
	protected AugmentedIntervalTreeNode insert(AugmentedIntervalTreeNode node, Interval s) {
		if(node != null) {
			int cmp = node.compare(s);
			// if it goes left, then s < node.str
			// this means that cmp > 0
			if(cmp <= 0) {
				node.left = this.insert((AugmentedIntervalTreeNode) node.left, s);
		//		node.update();
			} else if (cmp > 0) {
				node.right = this.insert((AugmentedIntervalTreeNode) node.right, s);
		//		node.update();
			}
			
			node = (AugmentedIntervalTreeNode) checkBalance(node, s);
		} else {
			node = new AugmentedIntervalTreeNode(s);
		}
		
		return node;
	}
	
	@Override
	public Set<Interval> findOverlaps(Interval interval) {
		Set<Interval> overlaps = new HashSet<Interval>();
		if(this.chromosome.equals(interval.getChromosome())) {	
			findOverlaps((AugmentedIntervalTreeNode) this.root, interval, overlaps);
		}
		return overlaps;
	}
	
	@Override
	public Set<Interval> findOverlaps(List<Interval> intervals) {
		Set<Interval> overlaps = new HashSet<Interval>();
		for(Interval i : intervals) {
			if(this.chromosome.equals(i.getChromosome())) {
				findOverlaps((AugmentedIntervalTreeNode) this.root, i, overlaps);
			}
		}		
		return overlaps;
	}
	
	protected void findOverlaps(AugmentedIntervalTreeNode node, Interval interval, Set<Interval> overlaps) {
		if(node != null) {			
			Interval contents = node.get();                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
			
			if(contents.overlaps(interval)) {
				overlaps.add(contents);
			}
			
			//findOverlaps((AugmentedIntervalTreeNode) node.getLeft(), interval, overlaps);
			//findOverlaps((AugmentedIntervalTreeNode) node.getRight(), interval, overlaps);
			
			if(node.getMax() > interval.getStart()) {
				findOverlaps((AugmentedIntervalTreeNode) node.getLeft(), interval, overlaps);
				
				if(contents.getStart() <= interval.getEnd()) {
					findOverlaps((AugmentedIntervalTreeNode) node.getRight(), interval, overlaps);
					if(contents.overlaps(interval)) {
						overlaps.add(contents);
					}
				}
			}
			
		}
	}
}
