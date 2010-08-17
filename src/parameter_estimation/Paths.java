package parameter_estimation;

import java.io.File;

/**
 * Paths contains paths to directories, files that are important in the Parameter Estimation program<BR>
 * Paths serves as a supertype to Param_Est, CKPackager, Rosenbrock types
 * @author nmvdewie
 *
 */
public class Paths {
	String workingDir;
	String chemkinDir;
	String binDir;
	String chemInp;
	String [] reactorInputs;
	String expDb;
	public String getExpDb() {
		return expDb;
	}
	String outputDir;

	//no_licenses sets the limiting number for the counting semaphore
	int noLicenses;
	
	public Paths (String workingDir, String chemkinDir, String chemInp, String [] reactorInputs, int noLicenses, String expDb){
		this.workingDir = workingDir;
		this.chemkinDir = chemkinDir;
		this.chemInp = chemInp;
		this.reactorInputs = reactorInputs;
		this.outputDir = workingDir+"output/";
		this.binDir = chemkinDir+"/bin/";
		/**
		 * TODO no. of licenses should not be part of Paths
		 */
		this.noLicenses = noLicenses;
		this.expDb = expDb;
		createOutputDir();
		
	}
	protected void createOutputDir (){
		boolean temp = new File(outputDir).mkdir();
		if(!temp){
			System.out.println("Creation of output directory failed!");
			System.exit(-1);
		}
	}

	public String getOutputDir() {
		return outputDir;
	}
	public void setOutputDir(String outputDir) {
		this.outputDir = outputDir;
	}
	public String getWorkingDir() {
		return workingDir;
	}
	public String getChemkinDir() {
		return chemkinDir;
	}
	public String getChemInp() {
		return chemInp;
	}
	public String[] getReactorInputs() {
		return reactorInputs;
	}
	public int getNoLicenses() {
		return noLicenses;
	}
	public String getBinDir() {
		return binDir;
	}
}
