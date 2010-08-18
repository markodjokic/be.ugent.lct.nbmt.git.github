package parameter_estimation;

import parameter_estimation.jing.mathtool.UncertainDouble;


/**
 * this class extends UncertainDouble (e.g. A,n,Ea) with a boolean 'fixed' that states whether this UncertainDouble
 * should be optimized (fixed = false) or not (fixed = true)
 * @author nmvdewie
 *
 */
public class FixedUncertainDouble {

	boolean fixed;
	
	public boolean isFixed() {
		return fixed;
	}
	public void setFixed(boolean fixed) {
		this.fixed = fixed;
	}
	UncertainDouble uDouble;
	
	public UncertainDouble getuDouble() {
		return uDouble;
	}
	public void setuDouble(UncertainDouble uDouble) {
		this.uDouble = uDouble;
	}
	/*
	 * standard value is false
	 */
	public FixedUncertainDouble(){
		this.fixed = false;
	}
	public FixedUncertainDouble(boolean fixed){
		this.fixed = fixed;
	}
		
}
