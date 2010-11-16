import jing.param.Temperature;
import parameter_estimation.ParameterEstimationDriver;




public class test {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String workingDir = "C:/Documents and Settings/nmvdewie/My Documents/Nick's Black Magic Toolbox/Eclipse Workspace/be.ugent.lct.nbmt.git.github/";
		String filename = "pinanol-raden_p1.inp";
		Temperature t = ParameterEstimationDriver.calcAvgTempPerExperiment(workingDir, filename);
		System.out.println(t.getK());
	}
		
		
}
