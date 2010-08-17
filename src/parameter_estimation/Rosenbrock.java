package parameter_estimation;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Map;
/**
 * Rosenbrock algorithm developed by Nick Vandewiele, February 2010 <BR>
 * inspired by B. Debrabandere's fortran implementation <BR>
 * adopted from "Computational Techniques for Chemical Engineers", by H.H. Rosenbrock and C. Storey <BR>
 * Ghent University <BR>
 * The Rosenbrock type performs the actual parameter optimization. Function type is called to calculate the Sum of Squared Errors (SSQ) <BR>
 * 
 * Rosenbrock algorithm uses two types as server classes:
 * 	-Optimization: provides getModelValues method (that calls update_chemistry_input to adjust chemistry input file) and 
 *   attributes like total_no_parameters, maxeval. Rosenbrock delegates this to Optimization type
 * 	-Function
 */

public class Rosenbrock{
	/**
	 * @param args
	 */
	public int maxeval;
	public double EFRAC = 0.3; //rosenbrock parameter
	public double SUCC = 3.0; //rosenbrock parameter
	public double FAIL = -0.5; //rosenbrock parameter
	public double[][] basis;
	public Double [] lower;
	public Double [] upper;

	public Double[] parms;

	public double[] e;
	
	/**
	 * Optimization type is used by Rosenbrock as server class, from which methods like Optimization.getModelValues are taken
	 * and from which attributes like all the 1D vectors are copied (throught the Optimization.getters)
	 * this solution is not that elegant as it requires to add (hard-code) Optimization attribute in new optimization routines.
	 * this requires to mess in the external source code...
	 * TODO: how to create an optimization code that does not require the addition of Optimization type in its source code? 
	 */
	private Optimization optim;
/**
 * constructor with standard values of EFRAC, SUCC, FAIL
 * @param o
 * @param maxeval
 */
	public Rosenbrock(Optimization o, int maxeval){
		this.optim = o;
		this.maxeval = maxeval;
		this.parms = new Double[Tools.retrieveFittedParameters(o.getCoefficients()).length];
		//put parameter to be fitted in Double array:
		System.arraycopy(Tools.retrieveFittedParameters(o.getCoefficients()), 0, parms, 0, parms.length);
	
		lower = Tools.retrieveLowerBounds(o.getCoefficients());
		upper = Tools.retrieveUpperBounds(o.getCoefficients());
		
		this.e = new double[parms.length];
		for (int i = 0; i < e.length; i++) {
				this.e[i] = Tools.retrieveFittedParameters(o.getCoefficients())[i] * EFRAC;	
		}
	}
	/**
	 * constructor in case adjustments to EFRAC, SUCC, FAIL are made
	 * @param o
	 * @param efrac
	 * @param succ
	 * @param fail
	 * @param maxeval
	 */
	public Rosenbrock(Optimization o, double efrac, double succ, double fail, int maxeval) {
		this(o,maxeval);
		this.EFRAC = efrac;
		this.SUCC = succ;
		this.FAIL = fail;
	}

	public Double[] returnOptimizedParameters() throws IOException, InterruptedException{
		//basis needs to be declared with the correct dimensions:
		basis = new double [parms.length][parms.length];

		PrintWriter out = new PrintWriter(new FileWriter("output_Rosenbrock.txt"));
		PrintWriter outSSQ = new PrintWriter (new FileWriter("SSQ_Rosenbrock.csv"));
		int neval = 1;
		out.println("Current evaluation no.: "+neval);

		// initialization of basis: OK
		basis = basisInit(basis);
		
		//evaluate model predictions with initial guesses:
		//flag_CKSolnList = true
		List<Map<String,Double>> model = optim.getModelValues(true);
		
		// function evaluation in initial point
		Function f = new Function(model,optim.getExp());
		
		//even in the initial point, one already has model values, error variance matrix can thus be taken, not just response variables
		double initial = f.getSRES();
		System.out.println("Initial SSQ: "+initial);
		double current = initial;
		System.out.println("Current value: "+current);
		outSSQ.println(neval+","+initial);
		
		// Set all flags to 2, i.e. no success has been achieved in direction i
		int [] flag = new int [parms.length];
		double [] d = new double [parms.length];
		flag = resetFlag(flag);
		d = resetD(d);
		
		Function fNew;
		// Main rosenbrock loop
		while (neval < optim.maxeval) {
			Double [] parms_trial = new Double[parms.length];
			
			for (int i = 0; i < parms.length; i++) {
					System.out.println("Beta: ");
					print(Tools.retrieveFittedParameters(optim.getCoefficients()));
					
					//Parameters are slightly changed in the direction of basisvector 'i' to beta_new(j), j=1..np
					for (int j = 0; j < parms.length; j++) {

						//new parameter trials:
						parms_trial[j] = parms[j] + e[i]*basis[j][i];
						
						//checking lower/upper bounds, and replacing value if crossed:
						parms_trial[j] = (parms_trial[j] < lower[j]) ? lower[j] : parms_trial[j];
						parms_trial[j] = (parms_trial[j] > upper[j]) ? upper[j] : parms_trial[j];	
					}
					
					System.out.println("Beta new (to be tested): ");
					print(parms_trial);
					
					//print new parameter guesses to file: 
					for (Double dummy : parms_trial) {
						out.print(dummy+", ");
					}
					out.println();
					
					//Number of evaluations is updated:
					neval++;
					
					out.println("Current evaluation no.: "+neval);
					System.out.println("Evaluation no. "+neval);
					
					//do update of chemistry input with new parameter trials:
					Tools.update_chemistry_input(optim.getPaths(), parms_trial, optim.getCoefficients());
					
					//model predictions with new parameter guesses is called:
					//set flag_CKSolnList to false to prevent calling the CKSolnList creator once again:
					//flag_CKSolnList = false
					model = optim.getModelValues(false);
				
					//Evaluate (value 'trial') cost function with new parameter guesses [beta_new(j)]
					fNew = new Function(model,optim.getExp());
					double trial = fNew.getSRES();
					
					out.println("Trial SSQ: "+trial);
					if(trial < current){
						out.println("Woohoo! trial < current!");
						out.println("Old SSQ: "+current);
						out.println("New SSQ: "+trial);
						outSSQ.println(neval+","+trial);
						
				
						//put successful parameter values of parms_trial in parms and in optimization.coefficients:
						parms = parms_trial.clone();
						optim.setCoefficients(Tools.SetListWithVector(parms, optim.getCoefficients()));
						
						current = trial;
						for (int j = 0; j < d.length; j++) {
							d[j] += e[j];
						}
						e[i] *= SUCC;
						// If flag(i) .EQ. 1, at least one success has occurred along direction i
						if (flag[i] == 2) {
							flag[i] = 1;
						}
					}
					else {
						out.println("Damn. trial SSQ > current SSQ...");
						outSSQ.println(neval+","+trial);
						e[i] *= FAIL;
						//If flag(i) == 0, at least one success has been followed by a least one failure in direction i.
						if (flag[i] == 1){
							flag[i] = 0;
						}
						//Test for condition "success followed by failure" for each direction.  If flag2 == 0, the test is positive.
						int flag2 = 0;
						for (int j = 0; j < flag.length; j++) {
							flag2 += flag[j];
						}
						if (flag2 == 0){
							basis = setNewBasis(d,basis);
							basis = gramschmidt(basis);
							flag = resetFlag(flag);
							d = resetD(d);
						}
					}
				
				
				//If number of evaluations exceeds maxeval, jump out of 'for' loop and evaluate once more:
				if (neval >= maxeval){
					i = parms.length;
				}
			}
		}
		
		out.close();
		outSSQ.close();
				
		return parms;
	}
	/**
	 * Initialization of the basis, taking unit vectors in every direction of the parameters

	 * @param bbasis
	 */
	public double[][] basisInit (double[][] bbasis) {
	for (int i = 0; i < bbasis[0].length; i++){
		for (int j = 0; j < bbasis[0].length; j++) {
			bbasis[i][j] = 0.0;
		}
		bbasis[i][i] = 1.0;
	}
	return bbasis;
	}

	/**
	* If successes have been found in multiple directions, this means that a new basis should be chosen<BR>
	* Basis vectors should be chosen in the joint direction of the successes<BR>
	* This method sets a new basis, but the actual implementation is still a mystery to me<BR>
	 * @param dd
	 * @param bbasis
	*/
	public double[][] setNewBasis (double[] dd, double[][] bbasis){
		
	for (int i = 0; i < dd.length; i++) {
		for (int j = 0; j < dd.length; j++) {
			bbasis[j][i] = dd[i] * bbasis[j][i];
		}
		for (int k = i+1; k < dd.length; k++) {
			for (int j = 0; j < dd.length; j++) {
				bbasis[j][i] = bbasis[j][i] + dd[k] * bbasis[j][k];	
			}
		}
	}
	return bbasis;
	}

	/**
	* 
	* @param bbasis
	* @return An orthonormal basis is derived and returned from the matrix basis using the Gram-Schmidt algorithm
	*/
	public double[][] gramschmidt (double[][] bbasis) {
	//  gram schmidt orthonormalization
	
	for(int  i = 0; i < bbasis[0].length; i++){
		for(int k = 0; k < i-1; k++){
			double scal_prod = 0.0;
			for (int j = 0; j < bbasis[0].length; j++){
				scal_prod = scal_prod + bbasis[j][i] * bbasis[j][k];
			}
			for (int j = 0; j < bbasis[0].length; j++) {
				bbasis[j][i] = bbasis[j][i] - scal_prod * bbasis[j][k];
			}
		}
	// calculation of norms of every basis vector: 
		double norm = 0.0;
		for (int j = 0; j < basis[0].length; j++){
			norm = norm + basis[j][i] * basis[j][i];
		}
	// normalization of new bases:          
		for (int j = 0; j < basis[0].length; j++){
			basis[j][i] = basis[j][i] / Math.sqrt(norm);
		}
	}
	//  nit = nit+1;
	     
	return bbasis;
	}

	/**
	* 
	* @param fflag
	 * @return all flags are set to 2, meaning i.e. no success has been achieved in direction i
	*/
	public int [] resetFlag (int[] fflag){
		for (int i = 0; i < fflag.length; i++) {
			fflag[i] = 2;
		}
		return fflag;
	}

	public double [] resetD (double[] dd){
		for (int i = 0; i < dd.length; i++) {
			dd[i] = 0.0;
		}
		return dd;
	}

	
	
	public static void print(Double [] d){
		for (Double i : d) {
			System.out.print(i+" ");
		}
		System.out.println();
	}
}