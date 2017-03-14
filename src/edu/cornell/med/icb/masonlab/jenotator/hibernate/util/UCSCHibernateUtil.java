package edu.cornell.med.icb.masonlab.jenotator.hibernate.util;

import java.util.HashMap;
import java.util.Map;

import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;
import org.hibernate.service.ServiceRegistry;
import org.hibernate.service.ServiceRegistryBuilder;

public class UCSCHibernateUtil {
	private static final Map<String, SessionFactory> sessionFactories = new HashMap<String, SessionFactory>();

	private static void buildSessionFactory(String dbName) {
		Configuration configuration = new Configuration();
		configuration.configure("edu/cornell/med/icb/masonlab/jenotator/hibernate/ucsc/config/hibernate.cfg.xml");
		configuration.setProperty("hibernate.connection.url", "jdbc:mysql://genome-mysql.cse.ucsc.edu/" + dbName);
		
		SessionFactory sessionFactory = configuration.buildSessionFactory();
		sessionFactories.put(dbName, sessionFactory);
		
	   // SessionFactory sessionfactory = configuration.buildSessionFactory(serviceRegistry);
		//SessionFactory sessionfactory = configuration.buildSessionFactory();
/*
		ServiceRegistry serviceRegistry = new ServiceRegistryBuilder()
				.configure("edu/cornell/med/icb/masonlab/jenotator/config/hibernate/ucsc/hibernate.cfg.xml")
				.applySetting("hibernate.connection.url", "jdbc:mysql://genome-mysql.cse.ucsc.edu/" + dbName)
				.buildServiceRegistry();
		SessionFactory sessionFactory = configuration.buildSessionFactory(serviceRegistry);
		sessionFactories.put(dbName, sessionFactory);
*/
	}

	public static SessionFactory getSessionFactory(String dbName) {
		if(!sessionFactories.containsKey(dbName)) {
			buildSessionFactory(dbName);
		}
		
		return sessionFactories.get(dbName);
	}
}
