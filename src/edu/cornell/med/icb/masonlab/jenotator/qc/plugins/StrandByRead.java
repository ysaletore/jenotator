package edu.cornell.med.icb.masonlab.jenotator.qc.plugins;

import java.util.HashMap;
import java.util.Map;

import edu.cornell.med.icb.masonlab.jenotator.model.MutableInt;
import net.sf.samtools.SAMRecord;

public class StrandByRead implements QCPlugin {
	
	private Map<String, Map<Integer, MutableInt>> positive;
	private Map<String, Map<Integer, MutableInt>> negative;

	public StrandByRead() {
		positive = new HashMap<String, Map<Integer, MutableInt>>();
		negative = new HashMap<String, Map<Integer, MutableInt>>();
	}
	
	@Override
	public void process(SAMRecord samrecord, String flowcell, int lane) {
		Map<Integer, MutableInt> flowcell_map;
		
		if(samrecord.getReadUnmappedFlag()) {
			flowcell_map = negative.get(flowcell);
			
			if(flowcell_map == null) {
				flowcell_map = new HashMap<Integer, MutableInt>();
				negative.put(flowcell, flowcell_map);
			}
		} else {
			flowcell_map = positive.get(flowcell);
			
			if(flowcell_map == null) {
				flowcell_map = new HashMap<Integer, MutableInt>();
				positive.put(flowcell, flowcell_map);
			}
		}
		
		MutableInt lane_counter = flowcell_map.get(lane);
		if(lane_counter == null) {
			lane_counter = new MutableInt();
			flowcell_map.put(lane, lane_counter);
		} else {
			lane_counter.increment();
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
