package edu.cornell.med.icb.masonlab.jenotator.model.hibernate;

import java.util.List;

import org.hibernate.Query;
import org.hibernate.Session;

import edu.cornell.med.icb.masonlab.jenotator.hibernate.ucsc.model.RefGene;
import edu.cornell.med.icb.masonlab.jenotator.hibernate.util.UCSCHibernateUtil;

public class BasicHibernateTest {
	public static void main(String[] args) {
		Session session = UCSCHibernateUtil.getSessionFactory("mm9").openSession();
		Query query = null;
		
		try {
			query = session.createQuery("from RefGene");
			@SuppressWarnings("unchecked")
			List<RefGene> list = (List<RefGene>) query.list();
			for(RefGene info : list) {
				System.out.println(info);
			}
		} catch(Throwable t) {
			t.printStackTrace();
		}
	}
}
