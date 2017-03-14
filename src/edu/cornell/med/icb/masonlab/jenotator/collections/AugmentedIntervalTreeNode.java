package edu.cornell.med.icb.masonlab.jenotator.collections;

import edu.cornell.med.icb.masonlab.jenotator.model.interval.Interval;


public class AugmentedIntervalTreeNode extends AVLNode<Interval> implements Comparable<AugmentedIntervalTreeNode> {
	protected int max;
	
	public AugmentedIntervalTreeNode(Interval s) {
		super(s);
		this.update();
	}

	public int getMax() {
		return this.max;
	}
	
	@Override
	public int compareTo(AugmentedIntervalTreeNode o) {
		if(this.contents.getStart() == this.contents.getStart()) {
			return ((Integer) this.contents.getEnd()).compareTo(o.contents.getEnd());
		} else {
			return ((Integer) this.contents.getStart()).compareTo(o.contents.getStart());
		}
	}
	
	@Override
	public void update() {
		super.update();
		
		if(this.left == null && this.right == null) {
			this.max = this.contents.getEnd();
		} else if(this.left == null) {
			((AVLNode<Interval>)this.right).update();
			this.max = Math.max(this.contents.getEnd(), ((AugmentedIntervalTreeNode) this.right).getMax());
		} else if(this.right == null) {
			((AVLNode<Interval>)this.left).update();
			this.max = Math.max(this.contents.getEnd(), ((AugmentedIntervalTreeNode) this.left).getMax());
		} else {
			((AVLNode<Interval>)this.left).update();
			((AVLNode<Interval>)this.right).update();
			this.max = Math.max(this.contents.getEnd(), 
					Math.max(((AugmentedIntervalTreeNode) this.left).getMax(), 
							((AugmentedIntervalTreeNode) this.right).getMax()));
		}
	}
}
