package levenberg.multi;

import java.io.BufferedReader;
import java.io.FileReader;

/** LMdemo Levenberg Marquardt demonstrator
  *
  * This one file contains four classes:
  *     LMdemo        contains just main()
  *     RosenHost     defines Rosenbrock demo problem.
  *     LMhost        interfaces methods needed by LM.
  *     LM            is the Levenberg-Marquardt routine.
  *
  *  M.Lampton UCB SSL (c) 1996; Java edition 2007
  */
public class LMdemomultiD
{
    public static void main(String args[]) throws Exception
    {

    	int no_experiments = 10;
    	int no_response_vars = 2;
    	int no_expl_vars = 2;
    	Double[][] exp = new Double[no_experiments][no_response_vars];
    	Double[][] explanatory = new Double[no_experiments][no_expl_vars];

    	
    	BufferedReader in = new BufferedReader(new FileReader("exp-trial.csv"));
    	for (int i = 0; i < no_experiments; i++){
    		String[] dummy = in.readLine().split(",");
    		for (int j = 0; j < no_expl_vars; j++)
	    		explanatory[i][j] = Double.parseDouble(dummy[j]);
    		
    		for (int k = 0; k < no_response_vars; k++)
    			exp[i][k] = Double.parseDouble(dummy[no_expl_vars+k]);    		
    	}
    	in.close();
    	
    	SimpleModelmultiDHost h = new SimpleModelmultiDHost(explanatory, exp);
    	
    	
    }
}
