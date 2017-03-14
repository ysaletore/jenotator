package edu.cornell.med.icb.masonlab.jenotator.activity;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.hibernate.Query;
import org.hibernate.Session;

import edu.cornell.med.icb.masonlab.jenotator.hibernate.ucsc.model.ChromInfo;
import edu.cornell.med.icb.masonlab.jenotator.hibernate.util.UCSCHibernateUtil;
import edu.cornell.med.icb.masonlab.jenotator.io.input.TabbedFileReader;

public class GetChromosomeSizes {
	public static Map<String, Integer> getDB(String genome) {
		Map<String, Integer> sizes = new HashMap<String, Integer>();
		Session session = UCSCHibernateUtil.getSessionFactory(genome).openSession();
		Query query = null;
		query = session.createQuery("from ChromInfo");
		@SuppressWarnings("unchecked")
		List<ChromInfo> list = (List<ChromInfo>) query.list();
		for(ChromInfo info : list) {
			sizes.put(info.getChrom(), (int) info.getSize());
		}
		
		return sizes;
	}
	
	public static Map<String, Integer> get(String filename) throws IOException {
		Map<String, Integer> sizes = new HashMap<String, Integer>();
		TabbedFileReader reader = new TabbedFileReader(filename);
		
		while(reader.hasNext()) {
			String[] line = reader.next();
			sizes.put(line[0], (int) Integer.parseInt(line[1]));
		}
		
		return sizes;
	}
}
