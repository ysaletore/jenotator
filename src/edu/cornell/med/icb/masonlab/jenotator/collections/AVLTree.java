package edu.cornell.med.icb.masonlab.jenotator.collections;

import java.util.Iterator;

/**
 * @author Yogesh Saletore
 * 
 * Creates an AVLTree
 */
public class AVLTree<T extends Comparable<T>> implements BinarySearchTree<T>  {
	protected BSTNode<T> root;
	protected int size;
	
	/**
	 * CONSTRUCTOR: AVLTree
	 *
	 */
	public AVLTree() {
		this.root = null;
		this.size = 0;
	}
	
	/**
	 * find: 				finds string s in tree
	 * @param 	T	s	string to search for
	 * @return	AVLNode		node that contains string,
	 *						null if couldn't be found
	 */
	public BSTNode<T> find(T s) {
		return find(this.root, s);
	}
	
	/**
	 * find: 					finds string s in tree
	 * @param 	AVLNode	node	node to start at
	 * @param 	String	s		string to search for
	 * @return	AVLNode			node that contains string,
	 * 							null if couldn't be found
	 */
	private BSTNode<T> find(BSTNode<T> node, T s) {
		if(node != null) {
			int cmp = node.compare(s);
			if(node.get() == s) {
				return node;
			} else if(cmp <= 0) {
				return find(node.left, s);
			} else {
				return find(node.right, s);
			}
		} else {
			return null;
		}
	}
	
	/**
	 * insert: 				insert string s into tree
	 * @param 	String	s	string to insert
	 */
	@Override
	public void insert(T s) {
		this.root = insert(this.root, s);
		this.size++;
	}
	
	/** 
	 * insert: 					insert string s into tree
	 * @param 	AVLNode	node	node to start insertion at
	 * @param	String	s		string to insert
	 * @return	AVLNode			the new node
	 */
	protected BSTNode<T> insert(BSTNode<T> node, T s) {
		if(node != null) {
			int cmp = node.compare(s);
			// if it goes left, then s < node.str
			// this means that cmp > 0
			if(cmp <= 0) {
				node.left = this.insert((AVLNode<T>)node.left, s);
			} else if (cmp > 0) {
				node.right = this.insert((AVLNode<T>)node.right, s);
			}
			node = checkBalance(node, s);
		} else {
			node = new AVLNode<T>(s);
		}
		return node;
	}
	
	/**
	 * checkBalance: ensures AVL structural condition is met
	 * @param	AVLNode	node	node to check
	 * @param	String	s		recently inserted string
	 * @return	AVLNode			new node
	 */
	protected BSTNode<T> checkBalance(BSTNode<T> n, T s) {
		AVLNode<T> node = (AVLNode<T>) n;
		node.update();			// reset height for current node
		
		int heightL, heightR;		// load the heights, in case null child
		if(node.left == null)
			heightL = -1;
		else
			heightL = ((AVLNode<T>)node.left).getHeight();
		if(node.right == null)
			heightR = -1;
		else
			heightR = ((AVLNode<T>)node.right).getHeight();
		
		int cmp1 = node.compare(s);
		
		/** if they differ by more than one, must rotate*/
		if(Math.abs(heightL - heightR) > 1) {
			if(cmp1 <= 0) {				// new one inserted at left
				int cmp2 = node.left.compare(s);	
				if(cmp2 <= 0) {			// single left rotate case
					BSTNode<T> a, b, x, y, z;
					a = node;
					b = a.left;
					z = a.right;
					x = b.left;
					y = b.right;
					
					b.left = x;
					b.right = a;
					a.left = y;
					a.right = z;
					
					((AVLNode<T>)a).update();
					((AVLNode<T>)b).update();
					return b;
				} else if (cmp2 > 0) {	// double left rotate
					BSTNode<T> a, b, c, w, x, y, z;
					
					a = node;
					b = a.left;
					z =  a.right;
					w =  b.left;
					c =  b.right;
					x =  c.left;
					y =  c.right;
					
					c.left = b;
					c.right = a;
					b.left = w;
					b.right = x;
					a.left = y;
					a.right = z;
					
					((AVLNode<T>)a).update();
					((AVLNode<T>)b).update();
					((AVLNode<T>)c).update();
					return c;
				}
			} else if (cmp1 > 0) {  // right
				int cmp2 = node.right.compare(s);
				if(cmp2 <= 0) {			// double rotate right
					BSTNode<T> a, b, c, w, x, y, z;
					
					a = node;
					b =  a.right;
					z =  a.left;
					w =  b.right;
					c =  b.left;
					x =  c.right;
					y =  c.left;
					
					c.right = b;
					c.left = a;
					b.right = w;
					b.left = x;
					a.right = y;
					a.left = z;
					
					((AVLNode<T>)a).update();
					((AVLNode<T>)b).update();
					((AVLNode<T>)c).update();
					return c;
				} else if (cmp2 > 0) {	// single right rotate
					BSTNode<T> a, b, x, y, z;
					a = node;
					b =  a.right;
					z =  a.left;
					x =  b.right;
					y =  b.left;
					
					b.right = x;
					b.left = a;
					a.left = z;
					a.right = y;
					
					((AVLNode<T>)a).update();
					((AVLNode<T>)b).update();
					return b;
				}
			}
		}
		
		node.update();
		return node;
	}
	
	/**
	 * getSize:		return size of tree
	 * @return	int	
	 */
	public int getSize() {
		return this.size;
	}
	
	@Override
	public String toString() {
		return this.root.toString();
	}

	@Override
	public Iterator<T> iterator() {
		return new BinarySearchTreeIterator<T>(this.root);
	}
	
	@Override
	public Iterator<T> iterator(BinarySearchTreeIterator.TraversalType type) {
		return new BinarySearchTreeIterator<T>(this.root, type);
	}
}
