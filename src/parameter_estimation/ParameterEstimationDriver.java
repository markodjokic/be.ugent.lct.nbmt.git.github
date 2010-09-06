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


public class ParameterEstimationDriver {

	/**
	 * @param args
	 * @throws Exception 
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		long time = System.currentTimeMillis();

		//input file will be searched in working directory under the name INPUT.txt:
		BufferedReader in = new BufferedReader(new FileReader(System.getProperty("user.dir")+"/INPUT.txt"));


		in.readLine();
		String workingDir = in.readLine();
		if (workingDir.charAt(workingDir.length()-1)!='/'){
			System.out.println("Pathname to working directory needs to end with forward slash '/' !");
			System.exit(-1);
		}

		in.readLine();
		String chemkinDir = in.readLine();

		in.readLine();
		int noLicenses = Integer.parseInt(in.readLine());

		in.readLine();
		String chemInp = in.readLine();

		in.readLine();
		String experimentsDb = in.readLine();

		in.readLine();
		boolean flagReactorDb = Boolean.parseBoolean(in.readLine());
		if (flagReactorDb){
			if(!(new File("reactor_input_template.inp").exists())){
				System.out.println("reactor_input_template.inp was not found in the working directory!");
				System.exit(-1);
			}
		}

		in.readLine();
		String reactorSetupsDb = in.readLine();
		
		in.readLine();
		Integer reactorSetupType = Integer.parseInt(in.readLine());
		
		in.readLine();
		int noExperiments = Integer.parseInt(in.readLine());

		//number of parameters to be fitted:
		in.readLine();
		in.readLine();
		
		//optimization flags:
		in.readLine();
		boolean flagRosenbrock = Boolean.parseBoolean(in.readLine());
		in.readLine();
		boolean flagLM = Boolean.parseBoolean(in.readLine());

		//maximum number of Rosenbrock evaluations:
		in.readLine();
		int maxeval = Integer.parseInt(in.readLine());

		//parameter estimation mode:
		in.readLine();
		int mode = Integer.parseInt(in.readLine());

		//REACTOR INPUT FILE NAMES:
		in.readLine();

		String [] reactor_inputs = new String[noExperiments];
		if (flagReactorDb){
			in.readLine();
		}
		else {
			for (int i = 0; i < noExperiments; i++){
				reactor_inputs[i] = in.readLine(); 
			}

		}


		//OPTIMIZATION SECTION: modified Arrhenius parameters [A, n, Ea]: 1 = loose, parameter will be varied; 0 = fixed, parameter will be fixed
		in.readLine();

		//number of reactions that will be optimized:
		in.readLine();
		int noFittedReactions = Integer.parseInt(in.readLine());

		List<ModifiedArrheniusKinetics> kin = new ArrayList();
		for (int i = 0; i < noFittedReactions; i++){
			//comment line with "reaction i: "
			in.readLine();
			// string of 1,0,1:
			String [] s = in.readLine().split(",");

			ConstrainedUncertainDouble A = new ConstrainedUncertainDouble(Integer.parseInt(s[0])==0);
			ConstrainedUncertainDouble n = new ConstrainedUncertainDouble(Integer.parseInt(s[1])==0);
			ConstrainedUncertainDouble Ea = new ConstrainedUncertainDouble(Integer.parseInt(s[2])==0);

			/**
			 * we populate the ArrayList kin, with the information we already have, i.e. information on whether or not
			 * some parameters need to be optimized.
			 */
			kin.add(new ModifiedArrheniusKinetics(A,n,Ea));
		}

		//		boolean flagNoParameters = checkNoParameters(fixRxns, noParams);
		boolean flagNoParameters = true;
		if (!flagNoParameters) {
			System.out.println("Number of parameters to be fitted provided in INPUT.txt does not equal the number of ones you specified in the optimization section in INPUT.txt!");
			System.exit(-1);
		}

		in.close();

		String template = "reactor_input_template.inp";

		/**
		 * for now, constraints on parameters are hard coded here:
		 */
		double min = 0;
		double max = 1e20;


		/**
		 * we populate min and max constraints in parameter container
		 */
		for(Iterator<ModifiedArrheniusKinetics> it = kin.iterator(); it.hasNext();){
			ModifiedArrheniusKinetics c = it.next();
			c.getA().setMin(min);
			c.getA().setMax(max);
			c.getN().setMin(min);
			c.getN().setMax(max);
			c.getEa().setMin(min);
			c.getEa().setMax(max);
		}
		
		if (flagReactorDb){
			if(reactorSetupType==1){
			reactor_inputs = reactorInputsParser(workingDir, reactorSetupsDb, template, noExperiments);
			}
			if(reactorSetupType==2){
				reactor_inputs = reactorInputsParser2(workingDir, reactorSetupsDb, template, noExperiments);
			}
		}
		
		
		Paths paths = new Paths(workingDir,chemkinDir,chemInp,reactor_inputs,noLicenses, experimentsDb);


		switch(mode){
		case 0:	System.out.println("PARITY PLOT MODE");	
		Param_Est p0 = new Param_Est(paths, kin, noExperiments, maxeval);
		p0.parity();
		Function f = new Function(p0.getModel(),p0.getExp());
		System.out.println("SSQ is: "+f.getSRES());
		break;
		
		case 1:	System.out.println("PARAMETER OPTIMIZATION MODE");		
		Param_Est p1 = new Param_Est(paths, kin, noExperiments, maxeval, flagRosenbrock, flagLM);
		p1.optimizeParameters();
		p1.statistics();
		p1.parity();
		break;
		
		case 2: System.out.println("EXCEL POSTPROCESSING MODE");
		Param_Est p2 = new Param_Est(paths);
		p2.excelFiles();
		break;
		
		case 3: System.out.println("STATISTICS MODE");
		Param_Est p3 = new Param_Est(paths, kin, noExperiments, maxeval, flagRosenbrock, flagLM);
		p3.statistics();
		p3.parity();
		}
		long timeTook = (System.currentTimeMillis() - time)/1000;
		System.out.println("Time needed for this program to finish: (sec) "+timeTook);
	}
	public static String[] reactorInputsParser(String workingDir, String experiments, String template, int no_experiments) throws IOException{
		ArrayList<String> reactorInputs = new ArrayList<String>();
		String filename;
		//read first line of excel input:
		/*
		 * Reading the first line will reveal:<BR>
		 * <LI>the number and names of species at the reactor inlet</LI>
		 * <LI>the T profile</LI>
		 * <LI>the P profile</LI>		
		 */
		BufferedReader in_excel = new BufferedReader(new FileReader(workingDir+experiments));
		String [] reactor_dim = in_excel.readLine().split(",");
		double convert_m_cm = 100;
		double length = Double.parseDouble(reactor_dim[1])*convert_m_cm;
		double diameter = Double.parseDouble(reactor_dim[2])*convert_m_cm;

		String [] dummy = in_excel.readLine().split(",");

		//		NOS : number of species
		int NOS = dummy.length-1;
		ArrayList<String> species_name = new ArrayList<String>();
		for (int i = 0; i < NOS; i++){
			species_name.add(dummy[1+i]);
		}

		dummy = in_excel.readLine().split(",");
		ArrayList<Double> species_mw = new ArrayList<Double>();
		for (int i = 0; i < NOS; i++){
			species_mw.add(Double.parseDouble(dummy[1+i]));
		}

		String [] exp = in_excel.readLine().split(",");

		int counter=0;
		while(!exp[counter].equals("TEMP")){
			counter++;
		}
		int position_T_profile = counter;
		counter++;

		//start reading the axial positions of the Temperature Profile:
		ArrayList<Double> array_Tprofile = new ArrayList<Double>();
		while(!exp[counter].equals("PRESS")){
			array_Tprofile.add(Double.parseDouble(exp[counter])*convert_m_cm);
			counter++;
		}
		int position_P_profile = counter;
		counter++;

		//start reading the axial positions of the Pressure Profile:
		ArrayList<Double> array_Pprofile = new ArrayList<Double>();
		for (int i = counter; i < exp.length; i++){
			array_Pprofile.add(Double.parseDouble(exp[i])*convert_m_cm);
		}
		/*
		 * Start writing the actual reactor input file, by doing the following:<BR>
		 * <LI>read in the lines of the template reactor input file that remain unchanged. write them to your output file</LI>
		 * <LI>change total mass flow rate</LI>
		 * <LI>add Pressure profile</LI>
		 * <LI>add Temperature profile</LI>
		 * <LI>add diameter</LI>
		 */

		int experiment_counter = 0;
		String line = null;
		try {
			line = in_excel.readLine();
			while(!dummy.equals(null)){				
				BufferedReader in_template = new BufferedReader(new FileReader(workingDir+template));
				String [] dummy_array = line.split(",");
				//experiment_counter contains the experiment number that will be used in the reactor input file name:
				experiment_counter = Integer.parseInt(dummy_array[0]);
				filename = "reactor_input_"+experiment_counter+".inp";
				
				reactorInputs.add(filename);
				PrintWriter out = new PrintWriter(new FileWriter(workingDir+filename));

				//copy the first 7 lines:
				for(int i = 0 ; i < 7 ; i++){
					String d = in_template.readLine();
					out.println(d);
				}

				//total mass flow rate:
				double massflrt = 0.0;
				double molflrt = 0.0;
				for(int i = 0 ; i < NOS; i++){
					massflrt=massflrt + Double.parseDouble(dummy_array[1+i])/3600;
					molflrt=molflrt + Double.parseDouble(dummy_array[1+i])/species_mw.get(i);
				}

				out.println("FLRT"+" "+massflrt);

				//Pressure Profile:
				double pressure = 0.0;
				double convert_bar_atm = 1.01325;
				for (int i = 0; i < array_Pprofile.size(); i++){
					pressure = Double.parseDouble(dummy_array[position_P_profile+i+1])/convert_bar_atm;
					out.println("PPRO "+array_Pprofile.get(i)+" "+pressure);
				}
				//Temperature Profile:
				double temperature = 0.0;
				double convert_C_K = 273.15;
				for (int i = 0; i < array_Tprofile.size(); i++){
					temperature = Double.parseDouble(dummy_array[position_T_profile+i+1])+convert_C_K;
					out.println("TPRO "+array_Tprofile.get(i)+" "+temperature);
				}

				//Diameter:
				out.println("DIAM "+diameter);

				//reactor length: 
				out.println("XEND "+length);
				//Inlet Species:
				double molfr = 0.0;
				for(int i = 0 ; i < NOS; i++){
					molfr = (Double.parseDouble(dummy_array[1+i])/species_mw.get(i))/molflrt;
					out.println("REAC "+species_name.get(i)+" "+molfr);
				}

				//force solver to use nonnegative species fractions:
				out.println("NNEG");

				//END:
				out.println("END");

				in_template.close();
				out.close();
				line = in_excel.readLine();
			}
			in_excel.close();
		}catch (Exception e){}//do nothing: e catches the end of the file exception

		// verify the correct number of lines in reactor input file:
		if( reactorInputs.size()!= no_experiments){
			System.out.println("Number of experiments in reactor inputs file does not correspond to the number of experiments provided in the INPUT file! Maybe check if .csv file contains redundant 'comma' lines.");
			System.exit(-1);
		}

		//convert ArrayList to String []:
		String[] a = new String[reactorInputs.size()];
		reactorInputs.toArray(a);
/*		String [] a = new String [reactorInputs.size()];
		for (int i = 0; i < reactorInputs.size(); i++){
			a[i] = reactorInputs.get(i);
		}
*/ 
		return a;
	}
	public static String[] reactorInputsParser2 (String workingDir, String experiments, String template, int no_experiments) throws IOException{
		ArrayList<String> reactorInputs = new ArrayList<String>();
		String filename;
		//read first line of excel input:
		/*
		 * Reading the first line will reveal:<BR>
		 * <LI>the number and names of species at the reactor inlet</LI>
		 * <LI>the T profile</LI>
		 * <LI>the P profile</LI>		
		 */
		BufferedReader in_excel = new BufferedReader(new FileReader(workingDir+experiments));

		String [] dummy = in_excel.readLine().split(",");

		//		NOS : number of species
		int NOS = dummy.length-1;
		ArrayList<String> species_name = new ArrayList<String>();
		for (int i = 0; i < NOS; i++){
			species_name.add(dummy[1+i]);
		}

		dummy = in_excel.readLine().split(",");
		ArrayList<Double> species_mw = new ArrayList<Double>();
		for (int i = 0; i < NOS; i++){
			species_mw.add(Double.parseDouble(dummy[1+i]));
		}


		/*
		 * Start writing the actual reactor input file, by doing the following:<BR>
		 * <LI>read in the lines of the template reactor input file that remain unchanged. write them to your output file</LI>
		 * <LI>change total mass flow rate</LI>
		 * <LI>add Pressure profile</LI>
		 * <LI>add Temperature profile</LI>
		 * <LI>add diameter</LI>
		 */

		int experiment_counter = 0;
		String line = null;
		in_excel.readLine();
		try {
			line = in_excel.readLine();
			
			double convert_m_cm = 100;
			double convert_C_K = 273.15;
			double diameter, length, temperature, pressure;
			
			while(!dummy.equals(null)){				
				BufferedReader in_template = new BufferedReader(new FileReader(workingDir+template));
				String [] dummy_array = line.split(",");
				//experiment_counter contains the experiment number that will be used in the reactor input file name:
				experiment_counter = Integer.parseInt(dummy_array[0]);
				filename = "reactor_input_"+experiment_counter+".inp";
				
				reactorInputs.add(filename);
				PrintWriter out = new PrintWriter(new FileWriter(workingDir+filename));

				//copy the first 7 lines:
				for(int i = 0 ; i < 7 ; i++){
					String d = in_template.readLine();
					out.println(d);
				}

				//total mass flow rate:
				double massflrt = 0.0;
				double molflrt = 0.0;
				for(int i = 0 ; i < NOS; i++){
					massflrt=massflrt + Double.parseDouble(dummy_array[1+i])/3600;
					molflrt=molflrt + Double.parseDouble(dummy_array[1+i])/species_mw.get(i);
				}

				out.println("FLRT"+" "+massflrt);

				convert_m_cm = 100;
				//Diameter:
				diameter = Double.parseDouble(dummy_array[NOS+1]) * convert_m_cm;
				out.println("DIAM "+diameter);

				//reactor length: 
				length = Double.parseDouble(dummy_array[NOS+2]) * convert_m_cm;
				out.println("XEND "+length);
				
				//temperature:
				temperature = Double.parseDouble(dummy_array[NOS+3])+convert_C_K;
				out.println("TPRO 0.0"+" "+temperature);
				out.println("TPRO "+length+" "+temperature);
				
				//pressure:
				pressure = Double.parseDouble(dummy_array[NOS+4]);
				out.println("PPRO 0.0"+" "+pressure);
				out.println("PPRO "+length+" "+pressure);
				
				//Inlet Species:
				double molfr = 0.0;
				for(int i = 0 ; i < NOS; i++){
					molfr = (Double.parseDouble(dummy_array[1+i])/species_mw.get(i))/molflrt;
					out.println("REAC "+species_name.get(i)+" "+molfr);
				}

				//force solver to use nonnegative species fractions:
				out.println("NNEG");

				//END:
				out.println("END");

				in_template.close();
				out.close();
				line = in_excel.readLine();
			}
			in_excel.close();
		}catch (Exception e){}//do nothing: e catches the end of the file exception

		// verify the correct number of lines in reactor input file:
		if( reactorInputs.size()!= no_experiments){
			System.out.println("Number of experiments in reactor inputs file does not correspond to the number of experiments provided in the INPUT file! Maybe check if .csv file contains redundant 'comma' lines.");
			System.exit(-1);
		}

		//convert ArrayList to String []:
		String [] a = new String[reactorInputs.size()];
		reactorInputs.toArray(a);

		return a;
	}
	public static boolean checkNoParameters(){
		return true;
	}
	
	
	

}
