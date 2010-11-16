import cern.jet.stat.Probability;


public class testProbability extends Probability {
	public static void main (String args[]){
		double alpha = 0.05; 
		int n = 100;
		int p = 4;
		double test = Probability.studentTInverse(alpha, n-p);
		System.out.println(test);
		
	
	}
	
}
