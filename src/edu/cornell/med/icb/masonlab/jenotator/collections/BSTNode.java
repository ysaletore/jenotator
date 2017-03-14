package edu.cornell.med.icb.masonlab.jenotator.collections;

public class BSTNode<T extends Comparable<T>> {
	protected BSTNode<T> left, right;
	protected T contents;

	public BSTNode() {
	}

	public BSTNode(T s) {
		contents = s;
		left = right = null;
	}

	public T get() {
		return this.contents;
	}

	/**
	 * compare: compare this node's contents with given contents
	 * @param 	T		s		object to compare to
	 * @return	int			value of comparison
	 * 			return = 0	objects are equal
	 * 			return > 0	given object goes to right
	 * 			return < 0	given object goes to left
	 */
	public int compare(T s) {
		return -this.contents.compareTo(s);
	}
	
	public BSTNode<T> getLeft() {
		return this.left;
	}
	
	public BSTNode<T> getRight() {
		return this.right;
	}

	@Override
	public String toString () {
		String out = "("+ contents.toString() + " ";
		if (left != null)
			out += left.toString();
		else out += ".";
		out += " ";
		if (right != null)
			out +=  right.toString();
		else out += ".";
		out += ")";
		
		out = "<";
		if(left != null) {
			out += left.toString();
		}
		out += ",";
		out += this.contents.toString();
		out += ",";
		if(right != null) {
			out += right.toString();
		}
		out += ">";
		
		return out;
	}
}