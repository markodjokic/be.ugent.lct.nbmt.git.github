package parameter_estimation;

import parsers.ConfigurationInput;
import datatypes.ModelValue;

public abstract class AbstractCKPackager extends Loggable{
	
	AbstractCKEmulation [] simulations;
	
	ConfigurationInput config;

	public ModelValue[] modelValues;
	
	public abstract AbstractCKEmulation []  runAllSimulations();

	public ConfigurationInput getConfig() {
		return config;
	}

	public ModelValue[] getModelValues() {
		return modelValues;
	}
}
