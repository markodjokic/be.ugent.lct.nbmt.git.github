package levenberg.original;


/**
 * This class will serve as a try-out to implement a simple non-linear regression with model 
 *  
 *  y1 = a * exp (-b * x1)
 *   
 *  
 * 	params a,b
 * 	explanatory var: x1
 * 	response vars: y1
 * @author nmvdewie
 *
 */

public class SimpleModel1DHost implements LMhost {

	//--------constants---------------
    protected final double SIGMA  = 1e-6;
    protected final double DELTAP = 1e-6;
    protected final double BIGVAL = 9.876543E+210; 
    protected int NPTS=10, NPARMS=2; 

    //--------fields------------------
    
    private double resid[] = new double[NPTS]; //residual matrix Y - Y^ with rows as experiments
    private double jac[][] = new double[NPTS][NPARMS];  
    private double parms[] = {1000.,10.};  // starting point
    
    private Double [] exp; //experimental data of response variables Y
    private Double [] model; //response matrix Y^
    private Double [] expl;//explanatory X matrix

    public SimpleModel1DHost(Double []x, Double[]exp) throws Exception
    {
    	this.exp = exp;
    	this.expl = x;
    	this.resid = new double[NPTS];
    	model = new Double[NPTS];
    	
    	for (int i=0; i<NPARMS; i++)
          System.out.println("Start parm["+i+"] = "+parms[i]); 

        LM myLM = new LM(this, NPARMS, NPTS); // run the minimizer

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
	            for (int i=0; i<NPTS; i++)
	              jac[i][j] = dGetResid(i);

	            for (int k=0; k<NPARMS; k++)
	              delta[k] = (k==j) ? -2*DELTAP : 0.0;

	            d = dNudge(delta); // resid at pminus
	            if (d==BIGVAL)
	            {
	                System.out.println("Bad dBuildJacobian() exit 3"); 
	                return false;  
	            }

	            for (int i=0; i<NPTS; i++)
	              jac[i][j] -= dGetResid(i);  // fetches resid[]

	            for (int i=0; i<NPTS; i++)
	              jac[i][j] *= FACTOR;

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
				resid[i] = (model[i]-exp[i]);
				sum += resid[i]*resid[i];
			 
		}
		
		return sum;
	}

	public double dGetJac(int i, int j) {
		 // Allows LM to get one element of the Jacobian matrix. 
	    
	        return jac[i][j]; 
	    }

	public double dGetResid(int i) {
		 // Allows LM to get one element of the resid[] vector. 
	    
	        return resid[i];
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
	public Double [] getModelValues(){
		for (int i = 0; i < model.length; i++){
			
    		model[i] = parms[0] * Math.exp(-parms[1] * expl[i]);  	
    		
    	}
		return model;
	}
}
