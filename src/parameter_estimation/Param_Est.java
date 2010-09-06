package parameter_estimation;
import java.util.*;
import java.io.*;

import jing.param.Temperature;




public class Param_Est{

	/**
	 * @throws InterruptedException 
	 */
	private Paths paths;
	public Paths getPaths() {
		return paths;
	}
	public int maxeval;
	
	public int noExp;
	private List<Map<String, Double>> exp;
	private List<Map<String,Double>> model;
	
	//private Parameters2D params;
	private List<ModifiedArrheniusKinetics> coefficients;
	public List<ModifiedArrheniusKinetics> getParams() {
		return coefficients;
	}
	public void setParams(List<ModifiedArrheniusKinetics> coefficients) {
		this.coefficients = coefficients;
	}
	
	//optimization flags:
	public boolean flagRosenbrock;
	public boolean flagLM;
	

	// constructor used for checking the validity of chemistry input file:
	public Param_Est(Paths paths){
		this.paths = paths;
	}
	//construct for parity mode:
	public Param_Est(Paths paths, List<ModifiedArrheniusKinetics> params, int no_experiments, int m_eval){
		this.paths = paths;
		this.coefficients = params;
		noExp = no_experiments;

		maxeval = m_eval;
	}
	
	//constructor used for parameter optimization option:
	public Param_Est(Paths paths, List<ModifiedArrheniusKinetics> params, int no_experiments, int m_eval, boolean flag_Rosenbrock, boolean flag_LM){
		this(paths, params, no_experiments, m_eval);
		this.flagRosenbrock = flag_Rosenbrock;
		this.flagLM = flag_LM;
	}
/**
 * optimizeParameters is the method that will optimize the kinetic parameters. It does so by:<BR>
 * <LI>checking the validity of the chemistry input file</LI>
 * <LI>taking the initial guesses of the kinetic parameters from the chemistry input file</LI>
 * <LI>reading the experimental database and store these values for the response variables</LI>
 * <LI>calling the actual optimization routine, i.e. the Rosenbrock algorithm</LI>
 * <LI>writing the optimized kinetic parameters to a params.txt file</LI>	
 * @throws Exception 
 * @throws Exception 
 */
	public void optimizeParameters() throws Exception{
		long time = System.currentTimeMillis();
		
		//check if initial input file is error-free:
		Runtime r = Runtime.getRuntime();
		CKEmulation c = new CKEmulation(paths, r);
		c.checkChemInput();
		
		// take initial guesses from chem.inp file:
		coefficients = Tools.initialGuess(paths.getWorkingDir(), paths.getChemInp(), coefficients);
		System.out.println("Initial Guesses of parameters are:");
		print(coefficients);
		
		//read experimental data:
		List<Map<String,Double>> exp = new ArrayList<Map<String,Double>>();
		exp = Tools.experimentsParser(paths.getExpDb(), noExp);
		//System.out.println(exp.toString());

		Optimization optimization = new Optimization(paths, coefficients, maxeval, flagRosenbrock, flagLM, exp);
		
		//call optimization routine:
		coefficients = optimization.optimize(exp);
		
		//print optimized parameters:
		System.out.println("New values of parameters are: ");
		print(coefficients);
		PrintWriter out = new PrintWriter(new FileWriter("params.txt"));

		for(ModifiedArrheniusKinetics con : coefficients){
			out.println("Reaction Parameters: ");
			out.print("A: "+con.getA().getuDouble().getValue()+"\\t");
			out.print("n: "+con.getN().getuDouble().getValue()+"\\t");
			out.println("Ea: "+con.getEa().getuDouble().getValue());
			out.println();
		}
		
		out.println();
		out.close();
		Tools.moveFile(paths.getOutputDir(), "params.txt");
		
		long timeTook = (System.currentTimeMillis() - time)/1000;
	    System.out.println("Time needed for this optimization to finish: (sec) "+timeTook);

	}
	/**	
 * plug new parameter guesses into the chemkin chemistry input file (usually chem.inp)
 * write new update chem.inp file
 * return the chemistry input filename
 * WARNING: method supposes a pre-conditioned chem.inp file, processed by ChemClean and with TD inside chem.inp!!!
 */
	public static String updateChemistryInput (String wd, double [] beta_new) throws IOException{
		String chemistry_input="chem.inp";
		String path_old_chem = wd+chemistry_input;	
		
		BufferedReader in = new BufferedReader(new FileReader(path_old_chem));
		PrintWriter out = new PrintWriter(new FileWriter(wd+"temp.inp"));
		String dummy = in.readLine();
		//just copy part of chem.inp about Elements, Species, Thermo
		boolean b = true;
		while(b){
			dummy = in.readLine();
			if (dummy.length() <= 8){
				b = true;
			}
			else if (dummy.substring(0,9).equals("REACTIONS")){
				b = false;
			}
			else {
				b = true;
			}
		}
		int i = 0;
		out.println(dummy);
		
		//read in next line, otherwise "REACTIONS" is still present in dummy String
		dummy = in.readLine();
		
		while(!dummy.equals("END")){
			String[] st_dummy = dummy.split("\\s");
			st_dummy[1] = Double.toString(beta_new[i]);
			i++;
			//System.out.println(i);
			st_dummy[3] = Double.toString(beta_new[i]);
			i++;
			//System.out.println(i);
			//toString(st_dummy);
			dummy = st_dummy[0];
			for (int j = 1; j < st_dummy.length; j++) {
				dummy = dummy +" "+st_dummy[j];
			}
			System.out.println(dummy);
			out.println(dummy);
			dummy = in.readLine();
		}
		// write "END" line:
		out.println(dummy);
		
		in.close();
		out.close();
		
		//delete old chem.inp file and create new chem.inp based on temp.inp:
		File f_old = new File(path_old_chem);
		f_old.delete();
		File f = new File(wd+"temp.inp");
		f.renameTo(new File(path_old_chem));
		return chemistry_input;
	}

	/**
	 * getModelPredictions_massfrac return a list with mass fractions as 'model values' instead of molar flowrates
	 * this becomes handy for parity plots, because experimental data will often be available in mass fractions, not molar flow rates
	 * @return
	 */
	public List<Map<String,Double>> getModelPredictionsMassfrac(){
		boolean flag_CKSolnList = true;
		CKPackager ckp = new CKPackager(paths, flag_CKSolnList);
		model = ckp.getModelValues();
		return model;
	}
	/**
	 * this routine produces model predictions without comparing them to experimental data
	 * @throws Exception 
	 */
	public void excelFiles() throws IOException, InterruptedException{
		long time = System.currentTimeMillis();
		//check if initial input file is error-free:
		Runtime r = Runtime.getRuntime();
		CKEmulation c = new CKEmulation(paths, r);
		c.checkChemInput();
		c.join();
			
		boolean flag_CKSolnList = true;
		boolean flag_toExcel = true;
		CKPackager ckp = new CKPackager(paths, flag_CKSolnList, flag_toExcel);
		ckp.getModelValues();
		
		moveOutputFiles();
		long timeTook = (System.currentTimeMillis() - time)/1000;
	    System.out.println("Time needed for Excel Postprocessing mode to finish: (sec) "+timeTook);
	}

	public void parity() throws IOException, InterruptedException{
		long time = System.currentTimeMillis();
		
		//check if initial input file is error-free:
		Runtime r = Runtime.getRuntime();
		CKEmulation c = new CKEmulation(paths, r);
		c.checkChemInput();
		c.join();
		
		List<Map<String,Double>> model = getModelPredictionsMassfrac();
		List<Map<String,Double>> exp = Tools.experimentsParser(paths.getExpDb(), noExp);
		List<String> speciesNames = Tools.getSpeciesNames(c.getPaths().getWorkingDir(), c.getAsu());

		//WRITE PARITY FILE:
		PrintWriter out = new PrintWriter(new FileWriter(paths.getWorkingDir()+"parity.csv"));
		
		// loop through all species:
		for(int i=0;i<speciesNames.size();i++){
			out.println(speciesNames.get(i).toString());
			// loop through all experiments:
			for(int j=0;j<exp.size();j++){
				Double experiment_value = exp.get(j).get(speciesNames.get(i));
				Double model_value = model.get(j).get(speciesNames.get(i));
				//out.println(speciesNames.get(i));
				out.println("experiment no. "+j+","+experiment_value+","+model_value+","+experiment_value);
			}
			out.println();
		}
		out.close();
		moveOutputFiles();
		Tools.moveFile(paths.getOutputDir(), "parity.csv");
		long timeTook = (System.currentTimeMillis() - time)/1000;
	    System.out.println("Time needed for Parity Mode to finish: (sec) "+timeTook);
	}
	
	public void print(List<ModifiedArrheniusKinetics> kin){
		for(ModifiedArrheniusKinetics c : kin){
			System.out.println("Reaction Parameters: ");
			System.out.println("A: "+c.getA().getuDouble().getValue());
			System.out.println("n: "+c.getN().getuDouble().getValue());
			System.out.println("Ea: "+c.getEa().getuDouble().getValue());
			System.out.println();
		}
	}
	public List<Map<String, Double>> getExp() throws IOException {
		exp = Tools.experimentsParser(paths.getExpDb(), noExp);
		return exp;
	}
	public List<Map<String, Double>> getModel() {
		return model;
	}
	public void statistics() throws Exception{
		long time = System.currentTimeMillis();
		
		//check if initial input file is error-free:
		Runtime r = Runtime.getRuntime();
		CKEmulation c = new CKEmulation(paths, r);
		c.checkChemInput();
		
		// take initial guesses from chem.inp file:
		coefficients = Tools.initialGuess(paths.getWorkingDir(), paths.getChemInp(), coefficients);
		System.out.println("Initial Guesses of parameters are:");
		print(coefficients);
		
		//read experimental data:
		List<Map<String,Double>> exp = new ArrayList<Map<String,Double>>();
		exp = Tools.experimentsParser(paths.getExpDb(), noExp);
		//System.out.println(exp.toString());

		Optimization optimization = new Optimization(paths, coefficients, maxeval, flagRosenbrock, flagLM, exp);
		
		optimization.calcStatistics();
		//moveOutputFiles();
		long timeTook = (System.currentTimeMillis() - time)/1000;
	    System.out.println("Time needed for this optimization to finish: (sec) "+timeTook);	    	    
	}
	protected void moveOutputFiles (){
		Tools.moveFiles(paths.getWorkingDir(), paths.getOutputDir(), ".out");
		Tools.moveFiles(paths.getWorkingDir(), paths.getOutputDir(), ".asu");
		Tools.moveFiles(paths.getWorkingDir(), paths.getOutputDir(), ".input");
		Tools.moveFiles(paths.getWorkingDir(), paths.getOutputDir(), ".asc");
		Tools.moveFile(paths.getOutputDir(),"CKSolnList.txt");
		
	}
	public Temperature getAvgTemp() throws Exception{
		return Tools.calcAvgTemp(paths.getWorkingDir(),paths.getReactorInputs());
	}
}
