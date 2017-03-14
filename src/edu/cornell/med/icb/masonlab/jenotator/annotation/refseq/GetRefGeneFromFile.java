package edu.cornell.med.icb.masonlab.jenotator.annotation.refseq;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.cornell.med.icb.masonlab.jenotator.hibernate.ucsc.model.RefGene;
import edu.cornell.med.icb.masonlab.jenotator.io.input.TextFileReader;

public class GetRefGeneFromFile {
	public static List<RefGene> load(String filename, String type) throws IOException {
		TextFileReader reader = new TextFileReader(filename);
		List<RefGene> refGenes = new ArrayList<RefGene>();
		String line;
		while(reader.hasNext() && ((line = reader.next()) != null)) {
			if(line.startsWith("#")) {
				continue;
			}
			RefGene refGene = new RefGene(line);
			refGenes.add(refGene);
		}
		
		return refGenes;
	}
}