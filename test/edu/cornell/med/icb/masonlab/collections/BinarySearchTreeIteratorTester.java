package edu.cornell.med.icb.masonlab.collections;

import java.util.Iterator;

import edu.cornell.med.icb.masonlab.jenotator.collections.AVLTree;
import edu.cornell.med.icb.masonlab.jenotator.collections.BinarySearchTree;

public class BinarySearchTreeIteratorTester {
	public static void main(String [] args) {
		BinarySearchTree<Integer> avl = new AVLTree<Integer>();
		int[] integers = {2,5,7,10,12,15,20,5,5,5,5,10,15,20,20,25};
		
		for(int i : integers) {
			avl.insert(i);
		}
		
		System.out.println(avl);
		
		// print traversals
		Iterator<Integer> in_order = avl.iterator();
		while(in_order.hasNext()) {
			System.out.print(in_order.next() + ", ");
		}
		System.out.println();
		
//		Iterator<Integer> pre_order = avl.iterator(BinarySearchTreeIterator.TraversalType.PRE_ORDER);
//		for(int i = pre_order.next(); pre_order.hasNext(); i = pre_order.next()) {
//			System.out.print(i + ", ");
//		}
//		System.out.println();
	}
}
