package levenberg.test;

import levenberg.original.LM;
import levenberg.original.LMhost;

public class ParaboloidHost implements LMhost {
	
	//--------constants---------------
    protected final double SIGMA  = 1e-6;
    protected final double DELTAP = 1e-6;
    protected final double BIGVAL = 9.876543E+210; 
	protected int NPTS=2, NPARMS=2; 
	
	//--------fields------------------
    
    private double resid[] = new double[NPTS]; 
    private double jac[][] = new double[NPTS][NPARMS];  
    private double parms[] = {-10.0, +10.0};  // starting point
	
	public ParaboloidHost () throws Exception{
	        for (int i=0; i<NPARMS; i++)
	          System.out.println("Start parm["+i+"] = "+parms[i]); 

	        LM myLM = new LM(this, NPARMS, NPTS); // run the minimizer

	        for (int i=0; i<NPARMS; i++)
	          System.out.println("End parm["+i+"]   = "+parms[i]); 
	}
	public boolean bBuildJacobian() {
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

	public double dComputeResid() {
		// Evaluates residual vector for parms[].
	    // Returns sum-of-squares. 
	    
	       
	        resid[0] = parms[0]; 
	        resid[1] = parms[1];
	        return resid[0]*resid[0] + resid[1]*resid[1]; 
	    

	}

	public double dGetJac(int i, int j) {
		// Allows LM to get one element of the Jacobian matrix. 
	    
	    return jac[i][j]; 
	}

	public double dGetResid(int i) {
		 // Allows LM to get one element of the resid[] vector. 
	    
	        return resid[i];
	}

	public double dNudge(double[] dp) {
		// Allows LM to modify parms[] and reevaluate.
	    // Returns sum-of-squares for nudged params.
	    // This is the only place that parms[] are modified.
	    // If NADJ<NPARMS, this is the place for your LUT.
	    
	        for (int j=0; j<NPARMS; j++)
	          parms[j] += dp[j]; 
	        return dComputeResid(); 
	}

}
