package parameter_estimation;

public class ConstrainedUncertainDouble extends FixedUncertainDouble{

	double min;
	public double getMin() {
		return min;
	}

	public void setMin(double min) {
		this.min = min;
	}
	
	double max;
	public double getMax() {
		return max;
	}

	public void setMax(double max) {
		this.max = max;
	}


	
	public ConstrainedUncertainDouble(double min, double max) {
		super();
		this.min = min;
		this.max = max;
	}

	public ConstrainedUncertainDouble(boolean fixed, double min, double max) {
		super(fixed);
		this.min = min;
		this.max = max;
	}

	public ConstrainedUncertainDouble(boolean fixed) {
		super(fixed);
	}
	
}
