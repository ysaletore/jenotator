package edu.cornell.med.icb.masonlab.jenotator.collections;

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Stack;

public class BinarySearchTreeIterator<T extends Comparable<T>> implements Iterator<T> {
	private TraversalType type;
	private Stack<BSTNode<T>> stack;
	private BSTNode<T> node;
	
	public enum TraversalType {
		IN_ORDER,
		PRE_ORDER,
		POST_ORDER;
	}
	
	public BinarySearchTreeIterator(BSTNode<T> root) {
		this(root, TraversalType.IN_ORDER);
	}
	
	public BinarySearchTreeIterator(BSTNode<T> root, TraversalType type) {
		this.type = type;
		this.stack = new Stack<BSTNode<T>>();
		
		switch(this.type){
			case PRE_ORDER:
				this.stack.push(root);
				break;
			case POST_ORDER:
				break;
			case IN_ORDER:
			default:
				this.node = root;
				break;
		}
	}

	@Override
	public boolean hasNext() {
		return (!this.stack.isEmpty() || this.node != null);
	}

	@Override
	public T next() {
		 if (!hasNext()) throw new NoSuchElementException();
		 T value;
		 
		 switch(this.type) {
		 	case PRE_ORDER:
		 		node = stack.pop();
		 		value = node.get();
			 	if(node.getRight() != null) stack.push(node.getRight());
			 	if(node.getLeft() != null) stack.push(node.getLeft());
			 	if(stack.isEmpty()) {
			 		node = null;
			 	}
			 	break;
			 	
		 	case POST_ORDER:
		 		throw new UnsupportedOperationException("POST_ORDER Iterator has not been implemented.");
//			 	break;
		 		
		 	case IN_ORDER:
		 	default:
			 	while(node != null) {
			 		stack.push(node);
			 		node = node.getLeft();
			 	}
			 	node = stack.pop();
			 	value = node.get();
			 	node = node.getRight();
			 	break;
		 }
		 
		 return value;
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}

}
