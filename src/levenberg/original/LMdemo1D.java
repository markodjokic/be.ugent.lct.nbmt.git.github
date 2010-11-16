package levenberg.original;

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
public class LMdemo1D
{
    public static void main(String args[]) throws Exception
    {
        //RosenHost myHost = new RosenHost();
    	//ParaboloidHost myHost = new ParaboloidHost();
    	int no_experiments = 10;
    	//int no_response_vars = 1;
    	int no_expl_vars = 1;
    	//Double[][] exp = new Double[no_experiments][no_response_vars];
    	Double[] exp = new Double[no_experiments];
    	//Double[][] expl = new Double[no_experiments][no_expl_vars];
    	Double[] explanatory = new Double[no_experiments];
    	
    	BufferedReader in = new BufferedReader(new FileReader("exp-trial.csv"));
    	for (int i = 0; i < no_experiments; i++){
    		String[] dummy = in.readLine().split(",");
	    		explanatory[i] = Double.parseDouble(dummy[0]);
    			exp[i] = Double.parseDouble(dummy[1]);    		
    	}
    	in.close();
    	
    	SimpleModel1DHost h = new SimpleModel1DHost(explanatory, exp);
    	
    	
    }
}
