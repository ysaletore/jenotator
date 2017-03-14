package edu.cornell.med.icb.masonlab.jenotator.collections;

import java.util.Iterator;

public interface BinarySearchTree<T extends Comparable<T>> extends Iterable<T> {
	public BSTNode<T> find(T s);
	public void insert(T s);
	public int getSize();
	public Iterator<T> iterator(BinarySearchTreeIterator.TraversalType type);
}
