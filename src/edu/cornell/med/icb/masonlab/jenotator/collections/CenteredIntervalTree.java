package edu.cornell.med.icb.masonlab.jenotator.collections;

import edu.cornell.med.icb.masonlab.jenotator.model.interval.Interval;

public class CenteredIntervalTree extends AVLNode<Interval> implements Comparable<AugmentedIntervalTreeNode> {

	public CenteredIntervalTree(Interval s) {
		super(s);
		// TODO Auto-generated constructor stub
	}

	@Override
	public int compareTo(AugmentedIntervalTreeNode o) {
		// TODO Auto-generated method stub
		return 0;
	}

}
