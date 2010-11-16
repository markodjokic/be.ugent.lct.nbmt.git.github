package levenberg.multi;
/**
 * This class will serve as a try-out to implement a multiresponse non-linear regression with model 
 *  
 *  y1 = a * exp (-b * x1)
 *  y2 = a * ln(b * x2)
 *  
 * 	params a,b
 * 	explanatory vars: x1,x2
 * 	response vars: y1,y2
 * @author nmvdewie
 *
 */

public class SimpleModelmultiDHost implements LMhostmultiD {

	//--------constants---------------
    protected final double SIGMA  = 1e-6;
    protected final double DELTAP = 1e-6;
    protected final double BIGVAL = 9.876543E+210; 
    protected int NPTS=10, NPARMS=2, NRESP = 2, NEXPL = 2; 

    //--------fields------------------
    
    private double resid[][] = new double[NPTS][NRESP]; //residual matrix Y - Y^ with rows as experiments
    private double jac[][][] = new double[NPTS][NPARMS][NRESP];  
    private double parms[] = {100.,100.};  // starting point
    
    private Double [][] exp = new Double [NPTS][NRESP]; //experimental data of response variables Y
    private Double [][] model = new Double [NPTS][NRESP];; //response matrix Y^
    private Double [][] expl = new Double [NPTS][NEXPL];//explanatory X matrix

    public SimpleModelmultiDHost(Double [][] x, Double[][]exp) throws Exception
    {
    	this.exp = exp;
    	this.expl = x;
    	
    	for (int i=0; i<NPARMS; i++)
          System.out.println("Start parm["+i+"] = "+parms[i]); 

        LMmultiD myLM = new LMmultiD(this, NPARMS, NPTS, NRESP); // run the minimizer

        for (int i=0; i<NPARMS; i++)
          System.out.println("End parm["+i+"]   = "+parms[i]); 
    }

	public boolean bBuildJacobian() throws Exception {
		// Allows LM to compute a new Jacobian.
	    // Uses current parms[] and two-sided finite difference.
	    // If current parms[] is bad, returns false.  
	    
	        double delta[] = new double[NPARMS];
	        double FACTOR = 0.5 / DELTAP; 
	        double d=0; 

	        for (int j=0; j<NPARMS; j++)
	        {
	            for (int k=0; k<NPARMS; k++)
	              delta[k] = (k==j) ? DELTAP : 0.0;

	            d = dNudge(delta); // resid at pplus
	            if (d==BIGVAL)
	            {
	                System.out.println("Bad dBuildJacobian() exit 2"); 
	                return false;  
	            }
	            for (int i=0; i<NPTS; i++){
	            	for (int k = 0; k < NRESP; k++){
	            		jac[i][j][k] = dGetResid(i,k);	
	            	}
	            } 

	            for (int k=0; k<NPARMS; k++)
	              delta[k] = (k==j) ? -2*DELTAP : 0.0;

	            d = dNudge(delta); // resid at pminus
	            if (d==BIGVAL)
	            {
	                System.out.println("Bad dBuildJacobian() exit 3"); 
	                return false;  
	            }
	            
	            for (int i=0; i<NPTS; i++){
	            	for (int k = 0; k < NRESP; k++){
	            		jac[i][j][k] -= dGetResid(i,k);	// fetches resid[]
	            	}
	            }
	            for (int i=0; i<NPTS; i++){
	            	for (int k = 0; k < NRESP; k++){
	            		jac[i][j][k] *= FACTOR;
	            	}
	            	
	            }
	            for (int k=0; k<NPARMS; k++)
	              delta[k] = (k==j) ? DELTAP : 0.0;

	            d = dNudge(delta);  
	            if (d==BIGVAL)
	            {
	                System.out.println("Bad dBuildJacobian() exit 4"); 
	                return false;  
	            }
	        }
	        return true; 
	}
	public double dComputeResid() throws Exception {
		getModelValues();
		Double sum = 0.0;
		for (int i = 0; i < resid.length; i++){
			for (int j = 0; j < resid[0].length; j++){
				resid[i][j] = (model[i][j]-exp[i][j]);
				sum += resid[i][j]*resid[i][j];
			}
		}
		
		return sum;
	}

	public double dGetJac(int i, int j, int k) {
		 // Allows LM to get one element of the Jacobian matrix. 
	    
	        return jac[i][j][k]; 
	    }

	public double dGetResid(int i, int j) {
		 // Allows LM to get one element of the resid[] vector. 
	    
	        return resid[i][j];
	    }

	public double dNudge(double[] dp) throws Exception {
		// Allows LM to modify parms[] and reevaluate.
	    // Returns sum-of-squares for nudged params.
	    // This is the only place that parms[] are modified.
	    // If NADJ<NPARMS, this is the place for your LUT.
	    
	        for (int j=0; j<NPARMS; j++)
	          parms[j] += dp[j]; 
	        return dComputeResid(); 
	}
	public Double [][] getModelValues(){
		for (int i = 0; i < model.length; i++){
			model[i][0] = parms[0] * Math.exp(-parms[1] * expl[i][0]);
    		model[i][1] = parms[0] * Math.log(parms[1] * expl[i][1]); //Math.log is natural log (base e) 	
			
    	}
		return model;
	}
	/**
	 * Jacobian with forward finite difference, to limit the number of extra "getModelValues()" calls,
	 * 1 instead of 2 for the two-sided finite difference method
	 */
	public boolean bBuildJacobian_forward() throws Exception {
		// Allows LM to compute a new Jacobian.
	    // Uses current parms[] and one-sided forward finite difference.
		// df/dx = ( f(x+h) - f(x) ) / h
	    // If current parms[] is bad, returns false.  
	    
	        double delta[] = new double[NPARMS];
	        double FACTOR = 1 / DELTAP; 
	        double d=0; 

	        for (int j=0; j<NPARMS; j++)
	        {	            
	        	for (int i=0; i<NPTS; i++){
	        		for (int k = 0; k < NRESP; k++){
	        			  jac[i][j][k] = -dGetResid(i,k); //resid central point
	        		}
	        	}   
     	
	        	for (int k=0; k<NPARMS; k++)
	              delta[k] = (k==j) ? parms[k]*DELTAP : 0.0;

	            d = dNudge(delta); // resid at pplus
	            if (d==BIGVAL)
	            {
	                System.out.println("Bad dBuildJacobian() exit 2"); 
	                return false;  
	            }
	            for (int i=0; i<NPTS; i++){
	            	for (int k = 0; k < NRESP; k++){
	            		jac[i][j][k] += dGetResid(i,k);	
	            	}
	              
	            }
	            for (int i=0; i<NPTS; i++){
	            	for (int k = 0; k < NRESP; k++){
	            		jac[i][j][k] *= FACTOR;	
	            	}
	            }
	            for (int k=0; k<NPARMS; k++)
		              delta[k] = (k==j) ? parms[k]*DELTAP : 0.0;

		        d = dNudge(delta);  
		        if (d==BIGVAL)
		        	{
		            System.out.println("Bad dBuildJacobian() exit 4"); 
		            return false;  
		        	}
	        }
	        return true; 
	}

	public double[][][] dGetFullJac() {
		return jac;
	}
}
