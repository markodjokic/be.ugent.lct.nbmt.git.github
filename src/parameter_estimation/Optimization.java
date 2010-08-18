package parameter_estimation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import stat.Statistics;

import levenberg.mono.NBMTHost;

/**Optimization type is to be seen as the 'driver' class of the broad set of optimization algorithms available, e.g. Rosenbrock.
 * Optimization type will condition the variables of NBMT in such a way that the actual optimization routine has to be modified to
 * serve for NBMT as little as possible. It does so, for instance, by converting the 2D parameter containers (double [][]) to 
 * 1D vectors (double[]), a format that is more natural to traditional optimization routines.
 * 
 * Also, the Optimization type will serve as a server class to the specific optimization routine (e.g. Rosenbrock) to provide methods
 * that access model values getModelValues, attributes e.g. total no_parameters, maximum number of evaluations
 * 
 * @author nmvdewie
 *
 */
public class Optimization{
	private Paths paths;
	public Paths getPaths() {
		return paths;
	}

	private List<ModifiedArrheniusKinetics> coefficients;

	public void setCoefficients(List<ModifiedArrheniusKinetics> coefficients) {
		this.coefficients = coefficients;
	}


	public List<ModifiedArrheniusKinetics> getCoefficients() {
		return coefficients;
	}

	public int maxeval;

	boolean flagRosenbrock;
	boolean flagLM;
	boolean weightedRegression;


	private List<Map<String,Double>> exp;

	/**
	 * Rosenbrock serves as server class to Optimization.
	 * it will execute the actual optimization
	 */
	private Rosenbrock rosenbrock;

	/**
	 * TODO name NBMTHost is not chosen very well... modify it! 
	 */

	private NBMTHost nbmthost;

	//constructor:
	public Optimization (Paths paths, List<ModifiedArrheniusKinetics> coefficients, int m, boolean f_r, boolean f_L, List<Map<String,Double>>exp){
		this.paths = paths;
		this.coefficients = coefficients;
		maxeval = m;

		flagRosenbrock = f_r;
		flagLM = f_L;
		//		weighted_regression = w_r;

		this.exp = exp;
	}


	public List<ModifiedArrheniusKinetics> optimize(List<Map<String,Double>> exp) throws IOException, InterruptedException{
		Set<String> response_vars = exp.get(0).keySet();
		PrintWriter out_species = new PrintWriter(new FileWriter("response_vars.txt"));
		for(Iterator<String> it = response_vars.iterator(); it.hasNext();){
			out_species.println((String)it.next());
		}
		out_species.close();

		if(flagRosenbrock){
/*			//Rosenbrock parameters:
			double efrac = 0.3;
			double succ = 3.0;
			double fail = -0.5;
*/			
			System.out.println("Start of Rosenbrock!");
			rosenbrock = new Rosenbrock(this, maxeval);
			rosenbrock.optimize();
			setCoefficients(Tools.setListWithVector(rosenbrock.getParms(), coefficients));
			
		}

		if(flagLM){
			System.out.println("Start of Levenberg-Marquardt!");
			nbmthost = new NBMTHost(this);
			setCoefficients(Tools.setListWithVector(nbmthost.getParms(), coefficients));		
		}


		// Rosenbrock monitors:
		if(new File("SSQ_Rosenbrock.csv").exists())
			Tools.moveFile(paths.getOutputDir(), "SSQ_Rosenbrock.csv");
		if(new File("output_Rosenbrock.txt").exists())
			Tools.moveFile(paths.getOutputDir(), "output_Rosenbrock.txt");

		//LM monitors:
		if(new File("LM.txt").exists())
			Tools.moveFile(paths.getOutputDir(), "LM.txt");
		if(new File("SSQ_LM.txt").exists())
			Tools.moveFile(paths.getOutputDir(), "SSQ_LM.txt");
		if(new File("response_vars.txt").exists())
			Tools.moveFile(paths.getOutputDir(), "response_vars.txt");

		return coefficients;
	}

/**
 * retrieves model values with the chemistry input file currently in place.
 * @param flag_CKSolnList
 * @return
 * @throws IOException
 * @throws InterruptedException
 */
	public List<Map<String,Double>> getModelValues(boolean flag_CKSolnList) throws IOException, InterruptedException{
		List<Map<String,Double>> model = new ArrayList<Map<String,Double>>();
		CKPackager ckp_new = new CKPackager(paths, flag_CKSolnList);
		model = ckp_new.getModelValues();
		return model;
	}

	public List<Map<String,Double>> getExp (){
		return exp;
	}
	/**
	 * converts the List<Map<String,Double>> format to a Double[][] format which is used in the LM optimization routine
	 * @return
	 */
	public Double[][] getExpDouble(){
		Double [][] dummy = new Double[exp.size()][exp.get(0).size()];

		// we want to have a fixed order in which the keys are called, therefore we put the response var names in a String []
		String [] species_names = new String [exp.get(0).size()];
		int counter = 0;
		for (String s : exp.get(0).keySet()){
			species_names[counter] = s;
			counter++;
		}

		for (int i = 0; i < exp.size(); i++){
			for (int j = 0; j < species_names.length; j++){
				dummy[i][j] = exp.get(i).get(species_names[j]);
			}
		}

		return dummy;
	}
	



	public NBMTHost getNBMTHost(){
		return nbmthost;
	}
	public void calcStatistics() throws IOException, InterruptedException{
		nbmthost = new NBMTHost(this, true);
		nbmthost.bBuildJacobian();
		Statistics s = new Statistics(this);
		PrintWriter out = new PrintWriter(new FileWriter("statistics.txt"));
		out.println("Averages of response variables:");
		out.println(this.nbmthost.getFunction().calcAverage());
		out.println();
		out.println("Variance-covariance of parameter estimations:");
		s.printMatrix(s.get_Var_Covar(), out);
		out.println();
		out.println("Correlation matrix of parameter estimations:");
		s.printMatrix(s.get_Corr(), out);
		out.println();
		out.println("t-values of individual significance of parameter estimations:");
		s.printArray(s.getT_values(), out);
		out.println();
		out.println("tabulated t-value for alpha = 5%");
		out.println(s.getTabulated_t_value());
		out.println();
		out.println("Confidence Intervals: [parameter][upper limit][lower limit]: ");
		s.printMatrix(s.getConfidence_intervals(), out);
		out.println();
		out.println("Number of experiments:");
		out.println(s.getNo_experiments());
		out.println("Number of parameters:");
		out.println(s.getNo_parameters());
		out.println("Number of responses:");
		out.println(s.getNo_responses());
		out.println("ANOVA: ");//Analysis of Variance:
		out.println("SRES: ");
		out.println(s.getSRES());
		out.println("SREG: ");
		out.println(s.getSREG());
		out.println("calculated F-value: ");
		out.println(s.getF_value());
		/*		
		out.println("tabulated F-value: ");
		out.println(s.getTabulated_F_value());
		 */		
		out.close();
		if(new File("statistics.txt").exists())
			Tools.moveFile(paths.getOutputDir(), "statistics.txt");
	}

	public boolean isWeighted_regression() {
		return weightedRegression;
	}

}
