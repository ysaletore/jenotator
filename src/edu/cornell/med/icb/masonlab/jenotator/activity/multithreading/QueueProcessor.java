package edu.cornell.med.icb.masonlab.jenotator.activity.multithreading;

import java.util.Queue;
import java.util.concurrent.Callable;

public class QueueProcessor<T> implements Callable<Integer> {
	protected final Queue<T> queue;
	protected final Processor<T> processor;
	protected final static int NUM_TRIES = 5;
	
	public QueueProcessor(Queue<T> queue, Processor<T> processor) {
		this.queue = queue;
		this.processor = processor;
	}

	@Override
	public Integer call() {
		try {
			Thread.sleep(1000);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		
		for(int i = 0; i < NUM_TRIES; i++) {
			while(!this.queue.isEmpty()) {
				T object = this.queue.remove();
				processor.process(object);
			}
			
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		return 1;
	}
	
	public Processor<T> getProcessor() {
		return this.processor;
	}
}
