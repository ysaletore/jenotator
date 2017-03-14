package edu.cornell.med.icb.masonlab.jenotator.io.input;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class TabbedFileReader {
	protected BufferedReader bufferedReader;
	protected String line;

	public TabbedFileReader(String filename) throws IOException {
		this(new File(filename));
	}
	
	public TabbedFileReader(File file) throws IOException {
		this.bufferedReader = new BufferedReader(new FileReader(file));
		this.line = this.bufferedReader.readLine();
	}
	
	public TabbedFileReader(InputStream inputstream) throws IOException {
		this.bufferedReader = new BufferedReader(new InputStreamReader(inputstream));
		this.readNext();
	}

	public String[] next() {
		if(this.line == null) {
			return null;
		}
		
		String[] parts = line.split("\t");
		try {
			this.readNext();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return parts;
	}

	public boolean hasNext() {
		return this.line != null;
	}
	
	private void readNext() throws IOException {
		this.line = this.bufferedReader.readLine();
	}
}
