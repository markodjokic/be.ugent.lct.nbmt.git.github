package parameter_estimation;

import java.io.*;
import java.util.*;
/**
 * CKEmulation is designed as a Thread type, implying that multiple CKEmulations can be initiated, allowing multithreading and possible speed-up<BR>
 * CKEmulation can call several Chemkin routines: Chem, CKReactorPlugFlow, GetSolution, GetSolnTranspose depending on the task to be executed<BR>
 * In order to cope with a limited number of Chemkin licenses, a counting Semaphore type is used to keep track of the number of licenses in use<BR>
 *  
 * @author nmvdewie
 *
 */
public class CKEmulation extends Thread{

	Paths paths;
	public Paths getPaths() {
		return paths;
	}

	private String reactorDir;
	
	private Runtime r;
	
	private String chem_out = "chem.out";
	private String chem_asc = "chem.asc";
	private String asu = "chem.asu";

	public String getAsu() {
		return asu;
	}



	// input file to CKPreProcess routine:
	private String preProc_inp = "CKPreProc_template.input";
	
	public String reactorSetup;
	public String reactorOut;
	
	public String CKSolnList = "CKSolnList.txt";
	public String xml = "XMLdata.zip";
	public String ckcsvName = "CKSoln.ckcsv";
	

	//species_fractions will be mol fractions or mass fractions depending on the flag massfrac. Anywho, the fractions are directly read through read_ckcsv
	public Map<String,Double> species_fractions;
	
	//'first' is flag that tells you if the CKSolnList needs to be constructed or not.
	boolean first;
	
	//'excel' tells you if the excel file (transposed CKSoln.ckcsv) needs to be created:
	boolean flagExcel;
	
	//Semaphore that controls chemkin license check-in and check-outs:
	Semaphore semaphore;
	
	//CONSTRUCTORS:
	//constructor for checking validity of chemistry input file:
	public CKEmulation(Paths paths, Runtime runtime){
		this.paths = paths;
		r = runtime;

	}
	
	//constructor for creating CKSolnList.txt
	public CKEmulation(Paths paths, Runtime runtime, String rs, boolean f){
		this(paths, runtime);
		reactorSetup = rs;
		first = f;
	}
	
	//constructor for running 'classical' Chemkin routines
	public CKEmulation(Paths paths, Runtime runtime, String rs, boolean f, Semaphore s, boolean excel){
		this(paths, runtime, rs, f);
		int length = rs.length();
		reactorOut = rs.substring(0,(length-4))+".out";
	
		this.flagExcel = excel;

		semaphore = s;
		
		this.reactorDir = paths.getWorkingDir()+"temp_ "+rs.substring(0,(length-4))+"/";
		boolean temp = new File(reactorDir).mkdir();
		if(!temp){
			System.out.println("Creation of reactor directory failed!");
			System.exit(-1);
		}
		
	}
	/**
	 * run() is the method that will be executed when Thread.start() is executed. Its argument list is void (mandatory I think).
	 */
	public void run(){
		try {
			semaphore.acquire();
			
			//System.out.println("license acquired!"+reactorSetup);
			
			//copy chem.inp to the reactorDir:
			Tools.copyFile(paths.getWorkingDir()+paths.getChemInp(),reactorDir+paths.getChemInp());
			Tools.copyFile(paths.getWorkingDir()+reactorSetup,reactorDir+reactorSetup);
			Tools.copyFile(paths.getWorkingDir()+"chemkindata.dtd",reactorDir+"chemkindata.dtd");
			callChem();	
			//call_PreProcess();
			callReactor();
			
			//copy reactor diagnostics file to workingdir:
			Tools.copyFile(reactorDir+reactorOut,paths.getWorkingDir()+reactorOut);
			
			//boolean first: if first time: create and adapt CKSolnList.txt file
			if (first){
				createSolnList();
				setSolnList();
				//copy the newly created CKSolnList to the workingDir so that it can be picked up by the other CK_emulations:
				Tools.copyFile(reactorDir+CKSolnList,paths.getWorkingDir()+CKSolnList);
			}
			else {
				//copy the CKSolnList to the reactorDir
				Tools.copyFile(paths.getWorkingDir()+CKSolnList,reactorDir+CKSolnList);
			}
			callGetSol();
			
			// if flag_excel = false: retrieve species fractions from the CKSoln.ckcsv file and continue:
			if (!flagExcel){
				//String name_data_ckcsv = "data_ckcsv"+j;

				species_fractions = Tools.readCkcsv(reactorDir,ckcsvName);
							
			}
			
			//if flag_excel = true: the postprocessed CKSoln.ckcsv file needs to be written to the parent directory (working directory)
			if (flagExcel){
				File excel_file = new File(reactorDir+ckcsvName);
				File dummy = new File (paths.getOutputDir()+ckcsvName+"_"+reactorSetup+".csv");
				excel_file.renameTo(dummy);
			}
			//delete complete reactorDir folder:
			Tools.deleteDir(new File(reactorDir));
			
			//when all Chemkin routines are finished, release the semaphore:
 			semaphore.release();
 			
			//System.out.println("license released!"+reactorSetup);
			
		} catch(Exception exc){
			System.out.println("Exception happened in CKEmulation run() method! - here's what I know: ");
			exc.printStackTrace();
			System.exit(-1);
			}
	}
	/**
	 * calls the .bat file that sets environment variables for proper use of future Chemkin calls<BR>
	 * @throws Exception
	 */
	public void callBat () throws IOException, InterruptedException {
		String [] setup_environment = {paths.getBinDir()+"run_chemkin_env_setup.bat"};
		executeCKRoutine(setup_environment);
	}
	
	/**
	 * call the Chemkin preprocessor chem and produces the linking file (.asc)
	 * @throws Exception
	 */
	public void callChem () throws IOException, InterruptedException {
		String [] preprocess = {paths.getBinDir()+"chem","-i",reactorDir+paths.getChemInp(),"-o",reactorDir+chem_out,"-c",reactorDir+chem_asc};
		executeCKRoutine(preprocess);
	}
	public void callPreProcess() throws IOException, InterruptedException {
		setCKPreProcess_Input();
		String [] preprocess = {paths.getBinDir()+"CKPreProcess","-i",paths.getWorkingDir()+preProc_inp};
		executeCKRoutine(preprocess);
	}
	/**
	 * call a Chemkin reactor model (ic: CKReactorPlugFlow) executable	
	 * @throws Exception
	 */
	public void callReactor () throws IOException, InterruptedException {
		String [] name_input_PFR = {paths.getBinDir()+"CKReactorPlugFlow","-i",reactorDir+reactorSetup,"-o",reactorDir+reactorOut};
		executeCKRoutine(name_input_PFR, new File(reactorDir));
	}
	
	/**
	 * processes the xmldata.zip file using Chemkin GetSolution executable
	 * @param flag_massfrac if true: GetSolution prints mass fractions instead of mole fractions (used for Excel Postprocessing and Parity mode)
	 * @throws Exception
	 */
	public  void callGetSol () throws IOException, InterruptedException {
		//String abbrev_path = cd+"data/abbreviations.csv";
		
		String [] progGetSol = {paths.getBinDir()+"GetSolution","-nosen","-norop","-mass",reactorDir+xml};
		executeCKRoutine(progGetSol, new File (reactorDir));
		
		Tools.deleteFiles(reactorDir, ".zip");
	}
	/**
	 * calls the Chemkin CKSolnTranspose executable
	 * @throws Exception
	 */
	public void callTranspose () throws Exception {
		String [] progTranspose = {paths.getBinDir()+"CKSolnTranspose",reactorDir+ckcsvName};
		executeCKRoutine(progTranspose);

	}
	
	 /**
	  * checkChemInput does a preliminary check of the initial chemistry output file to verify if no errors are present.<BR>
	  * It calls the Chemkin preprocessor which produces the output file<BR>
	  * This output file is read, and the String  " NO ERRORS FOUND ON INPUT: " is sought.<BR>
	  * If this String is not present, System.exit(-1) is called<BR>
	  */
	 public void checkChemInput(){
		 try {
			//create CKPreprocess.input file with directions to chem_inp, etc
			setCKPreProcess_Input();
			String [] preprocess = {paths.getBinDir()+"CKPreProcess","-i",paths.getWorkingDir()+preProc_inp};
			executeCKRoutine(preprocess);
				
				//read the produced chem.out (path_output) file, and check if it contains error messages:
				BufferedReader in = new BufferedReader(new FileReader(paths.getWorkingDir()+chem_out));
				String dummy = null;
				boolean flag = true;
				try {
					while(flag){
						dummy = in.readLine();
						if (dummy.equals(" NO ERRORS FOUND ON INPUT: ")){
							flag = false;
						}
					}
					in.close();
					if(!flag){
						System.out.println("Initial chemistry input file contains no errors. Proceed to parameter estimation!");
					}
					
				} catch(Exception e){
					System.out.println("Initial chemistry input file contains errors. Revision required!");
					System.exit(-1);
				}
		 }catch(Exception exc){
					System.out.println("exception happened - here's what I know: ");
					exc.printStackTrace();
					System.exit(-1);
		 }
	 }
	 /**
	  * createSolnList creates the CKSolnList.txt file by calling the "GetSolution -listonly" routine<BR>
	  * @throws Exception
	  */
	 public void createSolnList()throws Exception{
		String [] progGetList = {paths.getBinDir()+"GetSolution","-listonly",reactorDir+xml};
		executeCKRoutine(progGetList, new File(reactorDir));	
	 }
	 /**
	  * Set the SolnList.txt to the desired format:<BR>
	  * <LI>lines with # are left untouched</LI>
	  * <LI>no information whatsoever for all variables, except species, MW, exit_mass_flow_rate</LI>
	  * <LI>no sensitivity info for species, MW, exit_mass_flow_rate</LI>
	  * <LI>set UNIT of distance to (cm)</LI>
	  * <LI>all species mole fractions are reported, also those with negative fractions: FILTER MIN</LI>
	  * @throws Exception
	  */
	 public void setSolnList()throws Exception{
		 BufferedReader in = new BufferedReader(new FileReader(reactorDir+CKSolnList));
		 String temp = "tempList.txt";
		 PrintWriter out = new PrintWriter(new FileWriter(reactorDir+temp));
		 try{
			 String dummy = null;
			 dummy = in.readLine();
			 List<String> speciesNames = Tools.getSpeciesNames(paths.getWorkingDir(),asu);
			 //if a comment line (starts with char '#') is read, just copy it to output file
			 while(!dummy.equals(null)){
				 //if a comment line (#) or a blank line is read, just copy and continue
				 if (dummy.trim().length()==0){
					 out.println(dummy);
					 dummy = in.readLine();
					 //System.out.println(dummy);
				 }
				 else if (dummy.charAt(0)=='#'||(dummy.trim().length()==0)) {
					 out.println(dummy);
					 dummy = in.readLine();
					 //System.out.println(dummy);
				 }
				 
				 else {
					 //separator are TWO spaces, not just one space!
					 String[] st_dummy = dummy.split("  ");
					 //only species variables and molecular weight variable are reported:
					 if (st_dummy[0].equals("VARIABLE")){
						 //check if the 2nd keyword matches "molecular weight":
						 if (st_dummy[1].equals("molecular_weight")){
							//no sensitivity info for molecular weight variable
							 st_dummy[4]="0";
						 }
						 //check if the 2nd keyword matches "exit_mass_flow_rate":
						 else if(st_dummy[1].equals("exit_mass_flow_rate")){
							 st_dummy[4]="0";
						 }
						 //check if 2nd keyword is one of the species names:
						 else if(speciesNames.contains(st_dummy[1])){
							 //no sensitivity info for species: set last number in the line to zero
							 st_dummy[4]="0";
						 }
						 //the rest of the variables are set to zero and will not be reported in the .ckcsv file
						 else {
							 //st_dummy[3] is standard equal to zero
							 st_dummy[2]="0";
							 st_dummy[4]="0";
							 
						 }
					 }
					 //set UNIT of Distance to m instead of cm:
					 else if(st_dummy[0].equals("UNIT")){
						 if (st_dummy[1].equals("Distance")){
							 st_dummy[2]="(cm)";
						 }
					 }
					 
					 //make sure even negative mole fractions are printed:
					 else if(st_dummy[0].equals("FILTER")){
						 if (st_dummy[1].equals("MIN")){
							 st_dummy[2]="-1.0";
						 }
					 }
					 
					 //concatenate String array back to its original form:
					 String dummy_out = st_dummy[0];
					 //add double spaces between Strings again:
					 for(int i=1;i<st_dummy.length;i++){
						 dummy_out += "  "+st_dummy[i];
					 }
					 
					 out.println(dummy_out);
					 dummy = in.readLine();
					 //System.out.println(dummy);
				 }
			 } 
		}catch (Exception e){//do nothing: e catches the end of the file exception
			
			}
		in.close();
		out.close();
		File old = new File(reactorDir+CKSolnList);
		old.delete();
		File f_temp = new File(reactorDir+temp);
		f_temp.renameTo(new File(reactorDir+CKSolnList));
	}
	 
	 public void executeCKRoutine (String [] CKCommand) throws IOException, InterruptedException{
		 String s = null;		
		 Process p = r.exec(CKCommand);
			
		 BufferedReader stdInput_p = new BufferedReader(new InputStreamReader(p.getInputStream()));
		 BufferedReader stdError_p = new BufferedReader(new InputStreamReader(p.getErrorStream()));
	    
		// read the output from the command
	        //System.out.println("Here is the standard output of the command:\n");
	        while ((s = stdInput_p.readLine()) != null) {
	           //System.out.println(s);
	        }
	        stdInput_p.close();
	   // read any errors from the attempted command
	        //System.out.println("Here is the standard error of the command (if any):\n");
	        while ((s = stdError_p.readLine()) != null) {
	            //System.out.println(s);
	        }
	        stdError_p.close();
	        
			p.waitFor();
	        p.destroy();
			//System.out.println("Setup finished");		
	}
	 /**
	  * this routine overloads the standard execute_CKRoutine with a specified working directory, different from the standard working directory
	  * @param CKCommand
	  * @param workingDirectory
	  * @throws IOException
	  * @throws InterruptedException
	  */
	 public void executeCKRoutine (String [] CKCommand, File workingDirectory) throws IOException, InterruptedException{
		 String s = null;
		 String [] environment = null;
		 Process p = r.exec(CKCommand, environment, workingDirectory);
			
		 BufferedReader stdInput_p = new BufferedReader(new InputStreamReader(p.getInputStream()));
		 BufferedReader stdError_p = new BufferedReader(new InputStreamReader(p.getErrorStream()));
	    
		// read the output from the command
	        //System.out.println("Here is the standard output of the command:\n");
	        while ((s = stdInput_p.readLine()) != null) {
	            //System.out.println(s);
	        }
	        stdInput_p.close();
	   // read any errors from the attempted command
	        //System.out.println("Here is the standard error of the command (if any):\n");
	        while ((s = stdError_p.readLine()) != null) {
	            //System.out.println(s);
	        }
	        stdError_p.close();
	        
			p.waitFor();
	        p.destroy();
			//System.out.println("Setup finished");		
	}

	 
	 public Map<String,Double> getModelValue(){
		 /**
		  * TODO: this needs to be settled in a better (less patchy) way!
		  */
			 return species_fractions;

	 }


	 
	 public void setCKPreProcess_Input()throws IOException, InterruptedException{
		 //in windows: user.dir needs to be followed by "\", in *nix by "/"... 
		 String osname = System.getProperty("os.name");
		 if (osname.equals("Linux")){
			 PrintWriter out = new PrintWriter(new FileWriter(paths.getWorkingDir()+preProc_inp));
			 out.println("IN_CHEM_INPUT="+System.getProperty("user.dir")+"/"+paths.getChemInp());
			 out.println("OUT_CHEM_OUTPUT="+System.getProperty("user.dir")+"/"+chem_out);
			 out.println("OUT_CHEM_ASC="+System.getProperty("user.dir")+"/"+chem_asc);
			 out.println("OUT_CHEM_SPECIES="+System.getProperty("user.dir")+"/"+asu);
			 out.close();
		 }
		 else {
			 PrintWriter out = new PrintWriter(new FileWriter(paths.getWorkingDir()+preProc_inp));
		 out.println("IN_CHEM_INPUT="+System.getProperty("user.dir")+"\\"+paths.getChemInp());
		 out.println("OUT_CHEM_OUTPUT="+System.getProperty("user.dir")+"\\"+chem_out);
		 out.println("OUT_CHEM_ASC="+System.getProperty("user.dir")+"\\"+chem_asc);
		 out.println("OUT_CHEM_SPECIES="+System.getProperty("user.dir")+"\\"+asu);
		 out.close();
		 }
	 }
}


	

