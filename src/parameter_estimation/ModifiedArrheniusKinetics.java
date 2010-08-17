package parameter_estimation;

public class ModifiedArrheniusKinetics {
	public ConstrainedUncertainDouble A;
	
	public ConstrainedUncertainDouble getA() {
		return A;
	}

	public void setA(ConstrainedUncertainDouble a) {
		A = a;
	}
	
	public ConstrainedUncertainDouble n;
	public ConstrainedUncertainDouble getN() {
		return n;
	}

	public void setN(ConstrainedUncertainDouble n) {
		this.n = n;
	}
	
	public ConstrainedUncertainDouble Ea;
	public ConstrainedUncertainDouble getEa() {
		return Ea;
	}

	public void setEa(ConstrainedUncertainDouble ea) {
		Ea = ea;
	}

	public ModifiedArrheniusKinetics (ConstrainedUncertainDouble A, ConstrainedUncertainDouble n,ConstrainedUncertainDouble Ea){
		this.A = A;
		this.n = n;
		this.Ea = Ea;
		
	}
}
