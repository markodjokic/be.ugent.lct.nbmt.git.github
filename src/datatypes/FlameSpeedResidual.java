package datatypes;


public class FlameSpeedResidual extends Residual {

	public double value; 
	
	public FlameSpeedResidual(ExperimentalValue experimentalValue,
			ModelValue modelValue) {
		super.experimentalValue = (FlameSpeedExperimentalValue)experimentalValue;
		super.modelValue = (FlameSpeedModelValue)modelValue;
		compute();
	}

	@Override
	public void compute() {
		value = ((FlameSpeedModelValue)modelValue).value - ((FlameSpeedExperimentalValue)experimentalValue).value; 

	}

	@Override
	public Double getSSQValue() {
		return Math.pow(value,2);
	}

	@Override
	public double getValue() {
		return value;
	}

}
