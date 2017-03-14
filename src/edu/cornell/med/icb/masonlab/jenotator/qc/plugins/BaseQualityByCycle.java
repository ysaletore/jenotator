package edu.cornell.med.icb.masonlab.jenotator.qc.plugins;

import java.util.HashMap;
import java.util.Map;

import edu.cornell.med.icb.masonlab.jenotator.model.MutableInt;

import net.sf.samtools.SAMRecord;

public class BaseQualityByCycle implements QCPlugin {
	public static int QUALITY_FLOOR = 33;
	
	private Map<String, Map<Integer, Map<Integer, Map<Integer, MutableInt>>>> mapped;
	private Map<String, Map<Integer, Map<Integer, Map<Integer, MutableInt>>>> unmapped;
	
	public BaseQualityByCycle() {
		mapped = new HashMap<String, Map<Integer, Map<Integer, Map<Integer, MutableInt>>>>();
		unmapped = new HashMap<String, Map<Integer, Map<Integer, Map<Integer, MutableInt>>>>();
	}

	@Override
	public void process(SAMRecord samrecord, String flowcell, int lane) {
		byte[] base_quals = samrecord.getBaseQualities();
		Map<Integer, Map<Integer, Map<Integer, MutableInt>>> flowcell_map;
		
		if(samrecord.getReadUnmappedFlag()) {
			flowcell_map = unmapped.get(flowcell);
			
			if(flowcell_map == null) {
				flowcell_map = new HashMap<Integer, Map<Integer, Map<Integer, MutableInt>>>();
				unmapped.put(flowcell, flowcell_map);
			}
		} else {
			flowcell_map = mapped.get(flowcell);
			
			if(flowcell_map == null) {
				flowcell_map = new HashMap<Integer, Map<Integer, Map<Integer, MutableInt>>>();
				mapped.put(flowcell, flowcell_map);
			}
		}
		
		Map<Integer, Map<Integer, MutableInt>> lane_map = flowcell_map.get(lane);
		if(lane_map == null) {
			lane_map = new HashMap<Integer, Map<Integer, MutableInt>>();
			flowcell_map.put(lane, lane_map);
			
			for(int i = 0; i < base_quals.length; i++) {
				lane_map.put(i, new HashMap<Integer, MutableInt>());
			}
		}
		
		for(int i = 0; i < base_quals.length; i++) {
			MutableInt count = lane_map.get(i).get(base_quals[i]);
			if(count == null) {
				lane_map.get(i).put(new Integer(base_quals[i]), new MutableInt());
			} else {
				count.increment();
			}
		}
	}

	@Override
	public void exportPNG(String filename) {
		// TODO Auto-generated method stub

	}

	@Override
	public void exportPDF(String filename) {
		// TODO Auto-generated method stub

	}

	@Override
	public void exportHTML(String filename) {
		// TODO Auto-generated method stub

	}

}
