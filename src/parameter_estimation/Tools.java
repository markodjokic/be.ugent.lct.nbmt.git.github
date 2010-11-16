package parameter_estimation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import jing.mathTool.UncertainDouble;
import jing.param.GasConstant;
import jing.param.Temperature;
import jing.rxn.ArrheniusKinetics;


/**
 * Tools provides static methods to move/delete one File or multiple files with the same extension
 * @author nmvdewie
 *
 */
public class Tools {
	/**
	 * moveFiles moves all files with extension e to the specified destination dir
	 * @param orig_dir
	 * @param dest_dir
	 * @param e
	 */
	protected static void moveFiles( String orig_dir, String dest_dir, String e ) {
		ExtensionFilter filter = new ExtensionFilter(e);
		File original_dir = new File(orig_dir);
		String[] list = original_dir.list(filter);
		if (list.length == 0) return;
		boolean success;
		for (int i = 0; i < list.length; i++) {
			// Move file to new directory
			File file = new File(list[i]);
			success = file.renameTo(new File(dest_dir, list[i]));
			System.out.println( "File was successfully moved? " + success + "!");
		}    
	}

	protected static void moveFile(String dest_dir, String filename){
		// File (or directory) to be moved
		File file = new File(filename);

		// Move file to new directory
		boolean success = file.renameTo(new File(dest_dir, file.getName()));
	
	}
	/**
	 * from: http://forums.sun.com/thread.jspa?threadID=563148
	 * @param n
	 * @return
	 */
	public static void deleteDir(File dir){
		Stack<File> dirStack = new Stack<File>();
		dirStack.push(dir);

		boolean containsSubFolder;
		while(!dirStack.isEmpty()){
			File currDir = dirStack.peek();
			containsSubFolder = false;

			String[] fileArray = currDir.list();
			for(int i=0; i<fileArray.length; i++){
				String fileName = currDir.getAbsolutePath() + File.separator + fileArray[i];
				File file = new File(fileName);
				if(file.isDirectory()){
					dirStack.push(file);
					containsSubFolder = true;
				}else{
					file.delete(); //delete file
				}       
			}

			if(!containsSubFolder){
				dirStack.pop(); //remove curr dir from stack
				currDir.delete(); //delete curr dir
			}
		}
	}
	/**
	 * delete all files in directory d with extension e<BR>
	 * @param d
	 * @param e
	 */
	public static void deleteFiles( String d, String e ) {
		ExtensionFilter filter = new ExtensionFilter(e);
		File dir = new File(d);
		String[] list = dir.list(filter);
		File file;
		if (list.length == 0) return;

		for (int i = 0; i < list.length; i++) {
			file = new File(d + list[i]);
			boolean isdeleted =   file.delete();
			//System.out.print(file);
			//System.out.println( "  deleted " + isdeleted);
		}
	}

	public static List<ModifiedArrheniusKinetics> setListWithVector(Double [] vector, Optimization optim){
		//convert 1D vector back to List<ConstrainedArrheniusKinetics> notation:
		int counter = 0;
		List<ModifiedArrheniusKinetics> kin = optim.getCoefficients();
		/**
		 * first the List kin is updated with the values from Double[] vector, regardless of
		 * the fact whether re-parametrized A-values are used or not.
		 */
		for(ModifiedArrheniusKinetics c : kin){
			if(!c.getA().isFixed()) {
				//if(optim.isParametrized()){
					//how do we get the optimized value of Ea???
					//double A = vector[counter]*Math.exp(arg0);
					//c.getA().setuDouble(new UncertainDouble(vector[counter],0, "Adder"));
				//}
				c.getA().setuDouble(new UncertainDouble(vector[counter],0, "Adder"));
				counter++;
			}
			if(!c.getN().isFixed()){
				c.getN().setuDouble(new UncertainDouble(vector[counter],0, "Adder"));
				counter++;	
			}
			if(!c.getEa().isFixed()){
				c.getEa().setuDouble(new UncertainDouble(vector[counter],0, "Adder"));
				counter++;	
			}
		}

		/**
		 * if re-parametrization was used, the re-parametrized A-values need to be re-converted to
		 * the genuine A-values:
		 * A = k_avg * exp(Ea/(R * Tavg)) with k_avg stored in the A-container already.
		 */
		if(optim.isParametrized()){
			for(ModifiedArrheniusKinetics c : kin){
				if(!c.getA().isFixed()) {
					//since Ea is in kJ/mol, it has to be multiplied by 1000 to get J/mol:
					double Ea = c.getEa().getuDouble().getValue()*1000;
					double R = GasConstant.getJMolK();
					double Tavg = optim.getAvgTemp().getK();
					double kavg = c.getA().getuDouble().getValue(); 
					double A = kavg * Math.exp(Ea /(R * Tavg));
					c.getA().setuDouble(new UncertainDouble(A,0, "Adder"));
				}
			}
		}
		
		return kin;
	}
	/**
	 * convert_massfractions cuts of the Mass_fraction_ part off of the species names string
	 * @param m
	 * @return
	 */
	public static Map<String,Double> cutOffMassFrac_(Map<String,Double> m){
		Map<String, Double> dummy = new HashMap<String, Double> ();
		//loop through keys
		for ( String s : m.keySet()){
			//omit substring "Mass_fraction_" from key, i.e. take substring starting from character at position 14
			String dummy_name = s.substring(14);
			Double dummy_value = m.get(s);
			dummy.put(dummy_name, dummy_value);
		}
		return dummy;
	}	
	/**
	 * from: http://www.roseindia.net/java/beginners/CopyFile.shtml
	 * @param srFile
	 * @param dtFile
	 */
	public static void copyFile(String srFile, String dtFile){
		try{
			File f1 = new File(srFile);
			File f2 = new File(dtFile);
			InputStream in = new FileInputStream(f1);
			OutputStream out = new FileOutputStream(f2);

			byte[] buf = new byte[1024];
			int len;
			while ((len = in.read(buf)) > 0){
				out.write(buf, 0, len);
			}
			in.close();
			out.close();
			//System.out.println("File copied.");
		}
		catch(FileNotFoundException ex){
			System.out.println(ex.getMessage() + " in the specified directory.");
			System.exit(0);
		}
		catch(IOException e){
			System.out.println(e.getMessage());      
		}
	}
	/**
	 * read_ckcsv should read the CKSoln.ckcsv file and retrieve data from it.<BR>
	 * Which data specifically is explained here below:<BR>
	 * 	<LI>the mole fraction of all species</LI>
	 * the values should be taken at the end point of the reactor, i.e. the last data value of each row in the .ckcsv file<BR>
	 * the data will be stored in a LinkedList, chosen for its flexibility<BR>
	 * @param reactorDir directory where the file to be read is located
	 * @param ckcsv_name filename of file
	 * @throws IOException
	 */
	public static Map<String, Double> readCkcsv (String reactorDir, String ckcsv_name) throws IOException {
		Map <String, Double> dummy = new HashMap<String, Double>();
		BufferedReader in = new BufferedReader(new FileReader(reactorDir+ckcsv_name));

		String temp;
		String [] st_temp;
		LinkedList<String> list_temp;

		/*
		 *  Looking for the String "Exit_mass_flow_rate" since the line right after this line,
		 *  will contain the first species' mass fractions
		 */
		list_temp = new LinkedList<String>();
		do {
			list_temp.clear();
			temp = in.readLine();
			st_temp = temp.split(", ");
			for (int i=0;i<st_temp.length;i++){
				list_temp.add(st_temp[i]);
			}
		} while (!(list_temp.get(0)).equals("Exit_mass_flow_rate"));


		/* read all species' mass fractions, number of species is unknown, until the String "Molecular_weight" is encountered,
		 * which implies that the end of the list with species' mass fractions has been reached.
		 */
		list_temp.clear();
		list_temp = new LinkedList<String>();
		do {
			list_temp.clear();
			temp =in.readLine();
			st_temp = temp.split(", ");
			for (int i=0;i<st_temp.length;i++){
				list_temp.add(st_temp[i]);
			}
			if(!(list_temp.get(0)).equals("Molecular_weight")){				
				dummy.put((String)list_temp.get(0), Double.parseDouble(list_temp.get(list_temp.size()-1)));

			}
		} while (!(list_temp.get(0)).equals("Molecular_weight"));


		in.close();
		//convert to a HashMap with the real species names (cut of Mass_fraction_ or Mole_fraction_:
		dummy = cutOffMassFrac_(dummy);
		return dummy;
	}
	/**
	 * getSpeciesNames retrieves the names of the species from the chemistry input file
	 * @return
	 * @throws IOException
	 */
	public static List<String> getSpeciesNames(String workingDir, String chem_asu)throws IOException{
		BufferedReader in = new BufferedReader (new FileReader(workingDir+chem_asu));
		List<String> namesList = new ArrayList<String>();

		//first line contains number of species:
		int no_species = Integer.parseInt(in.readLine());

		String dummy = in.readLine();
		//first part of dummy contains: species=
		String dummy_speciesEq = dummy.substring(0, 8);
		//System.out.println(dummy_speciesEq);

		//rest of dummy contains species name and mw:
		String otherEnd = dummy.substring(8,dummy.length());
		//System.out.println(otherEnd);
		int index_mw = otherEnd.indexOf("mw=");
		//System.out.println(index_mw);
		String species_name = otherEnd.substring(0, index_mw-1).trim();
		//System.out.println(species_name);

		while(dummy_speciesEq.equals("species=")){
			namesList.add(species_name);
			dummy = in.readLine();
			if(dummy.length()>=8) {
				dummy_speciesEq = dummy.substring(0, 8);
				//System.out.println(dummy_speciesEq);

				//rest of dummy contains species name and mw:
				otherEnd = dummy.substring(8,dummy.length());
				//System.out.println(otherEnd);
				index_mw = otherEnd.indexOf("mw=");
				//System.out.println(index_mw);
				species_name = otherEnd.substring(0, index_mw-1).trim();
				//System.out.println(species_name);
			}
			else {
				break;
			}
		}


		in.close();
		if(no_species != namesList.size()){
			System.out.println("Something went wrong with the species names parsing from the .asu file!!!");
			System.exit(-1);
		}
		//System.out.println(namesList.toString());

		return namesList;
	}
	/**
	 * experiments_parser preprocesses the experimental data file into a easy-to-use format List<Map><String,Double>
	 * It reads the experimental data file. The format of this file (.csv) ought to be as follows:
	 * 	-first line states all the variable names, separated by commas
	 * 	-each of the following lines contain the experimental response variables, in the same order as the variable
	 *   names were declared
	 *  Each experiment is stored in a List with HashMaps as elements (see further)
	 *  The routines cuts each line in pieces, using the comma as separator
	 *  For each experiment, the variable name and the molar flow rate of each variable are put in a HashMap	 *   
	 * @return List of the experiments, with molar flowrates of the response variables
	 * @throws IOException
	 */
	public static List<Map<String, Double>> experimentsParser (String expName, int noExp)throws IOException{
		List<Map<String, Double>> exp = new ArrayList<Map<String, Double>>();
		try {
			//read experimental data file:
			BufferedReader in = new BufferedReader (new FileReader(expName));
			//read in species names on first line:
			String species_names = in.readLine();
			//System.out.println(species_names);
			String[] st_species = species_names.split(",");
			String dummy = in.readLine();
			//System.out.println(dummy);

			Map <String, Double> expMassFractions;
			while(dummy!=null){
				String[] st_dummy = dummy.split(",");
				expMassFractions = new HashMap <String, Double>();
				for (int j = 0; j < st_species.length; j++) {
					expMassFractions.put(st_species[j],Double.parseDouble(st_dummy[j]));	
				}
				exp.add(expMassFractions);
				//expMassFractions.clear();
				dummy = in.readLine();

			}
			if (exp.size()!= noExp){
				System.out.println("Experimental Database a different number of experiments than specified in INPUT file! Maybe check if .csv is created with redundand commas at the end...");
				System.exit(-1);
			}
			in.close();
		} catch(IOException e){
			System.out.println("Something went wrong during the preprocessing of the experimental data file!");
			System.exit(-1);
		}
		if((exp.size()==noExp)){
			return exp;
		}
		else{
			System.out.println("Experiments database contains different no. of experiments as defined in main class!");
			System.exit(-1);
			return null;	
		}

	}							
	/**
	 * initial_guess returns the initial parameter guesses, found in the chem.inp file.
	 * It does so by reading the file, searching the key-String "REACTIONS	KJOULES/MOLE	MOLES"
	 * from that point on, every line is read and the 2nd and 4th subString is taken and stored in a List l
	 * The 2nd and 4th element correspond to A and Ea of the modified Arrhenius equation
	 * The List l is then converted to a double array and returned
	 * @return initial guesses for kinetic parameters, as double array 
	 * @throws IOException
	 */
	public static List<ModifiedArrheniusKinetics> initialGuess (String workingDir, String chemInp, List<ModifiedArrheniusKinetics> kin) throws IOException{

		//			double[][] beta = new double[fixReactions.length][fixReactions[0].length];
		try {
			BufferedReader in = new BufferedReader(new FileReader(workingDir+chemInp));
			String dummy = in.readLine();

			//skip part of chem.inp about Elements, Species, Thermo
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

			/**
			 * 
			 * 			A new approach is taken, providing guidelines to the user to adapt his/her reaction mechanism file according to 
			 * 			what is specified in the guidelines:
			 * 			GUIDELINES:
			 * 			-use <=> to denote a reversible reaction (not =)
			 * 			-separate kinetic parameters by a single space
			 * 			-use a single space after the definition of the elementary reaction, e.g. A+B<=>C+D 1e13 0.0 200.0		
			 */
			for(ModifiedArrheniusKinetics c : kin){
				dummy = in.readLine();
				String[] st_dummy = dummy.split("\\s");
				//start with element at position 1 (not 0), because Arrhenius parameters start at position 1!
				c.getA().setuDouble(new UncertainDouble(Double.parseDouble(st_dummy[1]),0, "Adder"));
				c.getN().setuDouble(new UncertainDouble(Double.parseDouble(st_dummy[2]),0, "Adder"));
				c.getEa().setuDouble(new UncertainDouble(Double.parseDouble(st_dummy[3]),0, "Adder"));
			}

			in.close();

		} catch (IOException e) {
			System.out.println("Problem with obtaining initial parameter guesses!");
			System.exit(-1);
		}
		return kin;
	}
	/**
	 * stores the parameters to be optimized in a Double Array.
	 * It does so, by accessing the Coefficients
	 * first, it checks where a given parameter is fixed or not
	 * secondly, for pre-exponential factors, it checks whether the "parametrized" form
	 * is used for optimization instead of the "raw" value, which is usually very high
	 * e.g. 1E14 s-1 for certain types of unimolecular reactions 
	 * @param optim
	 * @return
	 */
	public static Double [] retrieveFittedParameters(Optimization optim){
		
		List<ModifiedArrheniusKinetics> l = optim.getCoefficients();
		List<Double> a = new ArrayList<Double>();
		for (ModifiedArrheniusKinetics m : l){			
			if(!m.getA().isFixed()){
				if(optim.isParametrized()){
					a.add(parametrize(m,optim.getAvgTemp()));
				}
				else a.add((Double)m.getA().getuDouble().getValue());
			}
			if(!m.getN().isFixed()) a.add((Double)m.getN().getuDouble().getValue());
			if(!m.getEa().isFixed()) a.add((Double)m.getEa().getuDouble().getValue());
		}

		Double [] dummy = a.toArray(new Double[a.size()]);
		return dummy;
	}

	public static Double [] retrieveLowerBounds(List<ModifiedArrheniusKinetics> l ){
		List<Double> a = new ArrayList<Double>();
		for (ModifiedArrheniusKinetics m : l){
			if(!m.getA().isFixed()) a.add((Double)m.getA().getMin());
			if(!m.getN().isFixed()) a.add((Double)m.getN().getMin());
			if(!m.getEa().isFixed()) a.add((Double)m.getEa().getMin());
		}
		Double [] dummy = a.toArray(new Double[a.size()]);
		return dummy;
	}
	
	public static Double [] retrieveUpperBounds(List<ModifiedArrheniusKinetics> l ){
		List<Double> a = new ArrayList<Double>();
		for (ModifiedArrheniusKinetics m : l){
			if(!m.getA().isFixed()) a.add((Double)m.getA().getMax());
			if(!m.getN().isFixed()) a.add((Double)m.getN().getMax());
			if(!m.getEa().isFixed()) a.add((Double)m.getEa().getMax());
		}
		Double [] dummy = a.toArray(new Double[a.size()]);
		return dummy;
	}
	
	/**	
	 * plug new parameter guesses into the chemkin chemistry input file (usually chem.inp)<BR>
	 * write new update chem.inp file <BR>
	 * return the chemistry input filename<BR>
	 * WARNING: method supposes TD inside chemistry input file!!!<BR>
	 */
	public static void writeChemistryInput (Paths p, List<ModifiedArrheniusKinetics> l) throws IOException{
		BufferedReader in = new BufferedReader(new FileReader(p.getWorkingDir()+p.getChemInp()));
		PrintWriter out = new PrintWriter(new FileWriter(p.getWorkingDir()+"temp.inp"));
		String dummy = in.readLine();

		//just copy part of chem.inp about Elements, Species, Thermo
		boolean b = true;
		while(b){
			out.println(dummy);
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
		out.println(dummy);
		
		//add new parameter values from Double[], only if isFixed()==true, otherwise, use old values from List
		for(ModifiedArrheniusKinetics m : l){
			dummy = in.readLine();
			String[] st_dummy = dummy.split("\\s");
			
			//put values of kinetic parameters (at position 1, 2, 3 of st_dummy[]):
			st_dummy[1] = Double.toString(m.getA().getuDouble().getValue());
			st_dummy[2] = Double.toString(m.getN().getuDouble().getValue());
			st_dummy[3] = Double.toString(m.getEa().getuDouble().getValue());
			
			//rebuild String by concatening, and write out:
			dummy = st_dummy[0];
			for (int j = 1; j < st_dummy.length; j++) {
				dummy = dummy +" "+st_dummy[j];
			}
			System.out.println(dummy);
			out.println(dummy);		
		}

		//just copy other reactions that are not varied, until end of file:
		dummy = in.readLine();
		while(!dummy.equals("END")){
			out.println(dummy);
			dummy = in.readLine();
		}

		out.println(dummy);

		in.close();
		out.close();

		File f_old = new File(p.getWorkingDir()+p.getChemInp());
		f_old.delete();
		File f = new File(p.getWorkingDir()+"temp.inp");
		f.renameTo(new File(p.getWorkingDir()+p.getChemInp()));
	}
	public static float calculateMachineEpsilonFloat() {
		float machEps = 1.0f;

		do {
			machEps /= 2.0f;
		}
		while ((float)(1.0 + (machEps/2.0)) != 1.0);

		return machEps;
	}
	/**
	 * reads reactor input file, searches for temperature profile, calculates average
	 * temperature along the profile
	 * @return
	 * @throws Exception 
	 */
	public static Temperature calcAvgTempPerExperiment (String workingDir, String filename) throws Exception{
		//open file:
		BufferedReader in = new BufferedReader(new FileReader(workingDir+filename));
		/**
		 * search for "TPRO" as first String in reactor input file line
		 */
		String line = in.readLine();
		while(!line.split(" ")[0].equals("TPRO")){
			line = in.readLine();
		}
		
		/**
		 * add all found temperatures to a dummy double; keep track of number of T points:
		 */
		double dummyTemp = 0.0;
		int counter = 0;
		while(line.split(" ")[0].equals("TPRO")){
			counter++;
			dummyTemp += Double.parseDouble(line.split(" ")[2]);
			line = in.readLine();
			
		}
		
		in.close();
		/**
		 * return average by dividing sum by number of T points: 
		 */
		return new Temperature(dummyTemp / counter,"K");
	}
	/**
	 * calculates average temperature:
	 * @param workingDir
	 * @param reactorInputs
	 * @return
	 * @throws Exception
	 */
	public static Temperature calcAvgTemp (String workingDir, String [] reactorInputs) throws Exception{
		double dummy = 0.0;
		for(String s : reactorInputs){
			dummy += Tools.calcAvgTempPerExperiment(workingDir, s).getK();
		}	
		return new Temperature (dummy / reactorInputs.length, "K");
	}
	/**
	 * convert pre-exponential factor A to reparametrized "average rate coefficient"
	 * @return
	 */
	public static Double parametrize(ModifiedArrheniusKinetics m, Temperature avgTemp){
		UncertainDouble A = m.getA().getuDouble();
		UncertainDouble n = m.getN().getuDouble();
		/**
		 * <code>ArrheniusKinetics</code> uses R in kcal/mol/K, therefore, we need to change the units of Ea to kcal/mol:
		 */
		double kJkCal = 0.239;
		UncertainDouble Ea = new UncertainDouble(m.getEa().getuDouble().getValue()*kJkCal,0, "Adder");
		ArrheniusKinetics arr = new ArrheniusKinetics(A,n,Ea,"Unknown",0,"Unknown","Unknown");
		return arr.calculateRate(avgTemp);
	}
	
	
}
