package parameter_estimation;

import datatypes.ModelValue;
import parsers.ConfigurationInput;
import readers.ReactorInput;

/**
 * Supertype for decorators for {@link CKEmulation}
 * @author nmvdewie
 *
 */
public abstract class CKEmulationDecorator extends AbstractCKEmulation {

	public AbstractCKEmulation simulation;
	
	@Override
	public abstract void run();

	public ConfigurationInput getConfig(){
		return simulation.getConfig();
	}
	

	public String getReactorDir() {
		return simulation.getReactorDir();
	}

	public ReactorInput getReactorInput() {
		return simulation.getReactorInput();
	}

	public String getReactorOut() {
		return simulation.getReactorOut();
	}
	
	public Runtime getRuntime() {
		return simulation.getRuntime();
	}
	
	public ModelValue getModelValue() {
		return simulation.getModelValue();
	}
	
}
