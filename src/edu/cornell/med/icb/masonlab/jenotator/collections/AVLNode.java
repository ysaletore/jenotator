package edu.cornell.med.icb.masonlab.jenotator.collections;

/**
 * Class: AVLNode
 * Internal class that creates a node for an AVLTree
 * @author Yogesh Saletore
 *
 */
public class AVLNode<T extends Comparable<T>> extends BSTNode<T> {
	/** height of node */
	protected int height;
	
	/**
	 * CONSTRUCTOR: AVLNode
	 * @param 		String	s		string to maintain count for
	 */
	public AVLNode(T s) {
		this.contents = s;
		this.left = this.right = null;
		this.height = 0;
	}
	
	/**
	 * resetHeight: resets the height for the current node
	 * 				Must be called after resetHeight is 
	 * 					called on children nodes.
	 */
	@SuppressWarnings("rawtypes")
	public void update() {
		if(this.left == null && this.right == null) {
			this.height = 0;
		} else if(this.left == null) {
			this.height = ((AVLNode)this.right).height + 1;
		} else if(this.right == null) {
			this.height = ((AVLNode)this.left).height + 1;
		} else {
			this.height = Math.max(((AVLNode)this.left).height, 
					((AVLNode)this.right).height) + 1;
		}
	}
	
	public int getHeight() {
		return this.height;
	}
	
	@Override
	public String toString() {
		return this.left + ", " + this.contents + ", " + this.right;
	}
}
