package parameter_estimation;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;

public class ExperimentalValues {
	private LinkedList<Map<String,Double>> experimentalEffluentValues;
	private LinkedList<Double> experimentalIgnitionValues;
	private LinkedList<Double> experimentalFlameSpeedValues;

	
	public ExperimentalValues(){
		experimentalEffluentValues = new LinkedList<Map<String,Double>>();
		experimentalIgnitionValues = new LinkedList<Double>();
	}

	/**
	 * calculates arithmetic average of all response variables, could be used as weights in regression
	 * @return
	 */
	public Map<String,Double> calcExperimentalEffluentAverage(){
		Map<String,Double> average = new HashMap<String,Double>();
		Set<String> response_vars = experimentalEffluentValues.get(0).keySet();
		for(Iterator<String> it = response_vars.iterator(); it.hasNext();){
			String s = (String) it.next();
			Double dummy = 0.0;
			for(int i = 0; i < experimentalEffluentValues.size(); i++){
				dummy = dummy + experimentalEffluentValues.get(i).get(s);
			}
			dummy = dummy / experimentalEffluentValues.size();
			average.put(s, dummy);
		}
		return average;
	}
	public Double calcExperimentalIgnitionAverage(){
		Double average = new Double(0.0);
		for(Iterator<Double> it = experimentalIgnitionValues.iterator(); it.hasNext();){
			average+=it.next();
		}
		return average;
	}
	public Double calcExperimentalFlameSpeedAverage(){
		Double average = new Double(0.0);
		for(Iterator<Double> it = experimentalFlameSpeedValues.iterator(); it.hasNext();){
			average+=it.next();
		}
		return average;
	}
	
	/**
	 * ####################
	 * GETTERS AND SETTERS:
	 * ####################
	 */
	
	/**
	 * @category getter
	 */
	public LinkedList<Double> getExperimentalFlameSpeedValues() {
		return experimentalFlameSpeedValues;
	}
	/**
	 * @category setter
	 */
	public void setExperimentalFlameSpeedValues(
			LinkedList<Double> experimentalFlameSpeedValues) {
		this.experimentalFlameSpeedValues = experimentalFlameSpeedValues;
	}
	/**
	 * @category getter
	 */
	public LinkedList<Map<String, Double>> getExperimentalEffluentValues() {
		return experimentalEffluentValues;
	}
	/**
	 * @category setter
	 */
	public void setExperimentalEffluentValues(
			LinkedList<Map<String, Double>> experimentalEffluentValues) {
		this.experimentalEffluentValues = experimentalEffluentValues;
	}
	/**
	 * @category getter
	 */
	public LinkedList<Double> getExperimentalIgnitionValues() {
		return experimentalIgnitionValues;
	}
	/**
	 * @category setter
	 */
	public void setExperimentalIgnitionValues(
			LinkedList<Double> experimentalIgnitionValues) {
		this.experimentalIgnitionValues = experimentalIgnitionValues;
	} 
}
