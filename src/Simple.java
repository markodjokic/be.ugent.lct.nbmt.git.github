// simple class to test the JUnit test framework
public class Simple {
	public static int sum (int[] vec){
		int total=0;
		for (int i = 0; i < vec.length; i++) {
			total += vec[i];	
		}
		return total;
	}
	public static int prod (int[]vec){
		int total=1;
		for (int i = 0; i < vec.length; i++) {
			total *= vec[i];	
		}
		return total;
	}
 
}
