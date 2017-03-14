package edu.cornell.med.icb.masonlab.jenotator.io.input;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class TextFileReader {
	protected BufferedReader bufferedReader;
	protected String line;

	public TextFileReader(String filename) throws IOException {
		this(new File(filename));
	}
	
	public TextFileReader(File file) throws IOException {
		this.bufferedReader = new BufferedReader(new FileReader(file));
		this.readNext();
	}
	
	public TextFileReader(InputStream inputstream) throws IOException {
		this.bufferedReader = new BufferedReader(new InputStreamReader(inputstream));
		this.readNext();
	}

	public String next() {
		if(this.line == null) {
			return null;
		}
		
		String line2 = line;
		try {
			this.readNext();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return line2;
	}

	public boolean hasNext() {
		return this.line != null;
	}
	
	private void readNext() throws IOException {
		this.line = this.bufferedReader.readLine();
	}
}
