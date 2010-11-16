
import cern.jet.stat.*;
public class testNickColt {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

/*		DoubleMatrix2D matrix;
		matrix = new DenseDoubleMatrix2D(3,4);
		System.out.println(matrix);
		
		int rows = matrix.rows();
		int cols = matrix.columns();
		
		double counter = 0;
		for (int i = 0; i < matrix.rows(); i++) {
			for (int j = 0; j < matrix.columns(); j++) {
				matrix.set(i, j, counter);
				counter++;
			}
		}
		
		System.out.println(matrix);
		
		//transpose:
		matrix = matrix.viewDice().copy();
		System.out.println(matrix);
		
		//we build a matrix we want to decompose:
		double[][] m = {{3,-1,-1},{-1,3,-1},{-1,-1,3}};
		matrix = new DenseDoubleMatrix2D(m);
		System.out.println(matrix);
		
		//Eigenvalue decomposition:
		EigenvalueDecomposition e = new EigenvalueDecomposition(matrix);
		DoubleMatrix2D D = e.getD();
		DoubleMatrix2D V = e.getV();
		System.out.println(D); //diagonal
		System.out.println(V);//eigenvectors
		
		//matrix multiplication occurs through the Algebra type:
		Algebra a = new Algebra();
		DoubleMatrix2D matrixbis = a.mult(a.mult(V,D),V.viewDice().copy());
		System.out.println("Product is:");
		System.out.println(matrixbis);
		
*/
	
		
	}

}
