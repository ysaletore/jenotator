package edu.cornell.med.icb.masonlab.collections;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class MultiCoreTester {
	public static class MyThread extends Thread implements Runnable {
		public void run() {
			while(true) {
				
			}
		}
	}
	
	public static void main(String [] args) {
		ExecutorService svc = Executors.newFixedThreadPool(40);
		for(int i = 0; i < 40; i++) {
			svc.submit(new MyThread());
		}
		
		while(true) {
		}
	}
}
