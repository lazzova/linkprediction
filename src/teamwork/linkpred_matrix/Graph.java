package teamwork.linkpred_matrix;

import java.util.Date;

import org.apache.commons.math3.random.GaussianRandomGenerator;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomVectorGenerator;
import org.apache.commons.math3.random.UncorrelatedRandomVectorGenerator;

public class Graph {
	 /** Theree dimensional matrix, each (i,j)th element is 
	  *  a vector of features
	  */
	
	/*
	public static double [] generateFeatures (int f) {
		/** Generate array of f Gaussian features *//*
		JDKRandomGenerator rand = new JDKRandomGenerator();
		rand.setSeed(new Date().getTime()); 
		RandomVectorGenerator randomVector = 
				new UncorrelatedRandomVectorGenerator(
				f, new GaussianRandomGenerator(rand));
		
		return randomVector.nextVector();
	}
	*/
	public static JDKRandomGenerator rand = new JDKRandomGenerator();
	public static RandomVectorGenerator randomVector = null;
	
	public static void initRandom (int f) {
		/** Set seed and initialize the generator of
		 *  random vectors of length f
		 */
		rand.setSeed(new Date().getTime());
		randomVector = 
				new UncorrelatedRandomVectorGenerator(
				f, new GaussianRandomGenerator(rand));
	}
	
	
	/*
	public static double [] generateFeatures () {
		/** Generate array of f Gaussian features *		
		return randomVector.nextVector();
	}
    */

	public static byte [][] generate (int n) {
		/** Generate undirected graph with n nodes.
		 *  Used for testing purpose		 
		 */
		
		//JDKRandomGenerator rand = new JDKRandomGenerator();
		//rand.setSeed(new Date().getTime()); 
				
		byte [][] g = new byte [n][n];		
		int [] degCumulative = new int [n];	                     // array for cumulative degree sums
						
		g[0][1] = 1;                                             // connect first three nodes in a triad
		g[1][0] = 1;
		g[1][2] = 1;  
		g[2][1] = 1; 
		g[2][0] = 1; 
		g[0][2] = 1; 
		
		int k;
		int len;
		int randNum;
		for (int i = 3; i < n; i++) {
			
			for (int j = 0; j < 3; j++) {                        // generate three links
				if (rand.nextInt(11) < 8) {                      // select destination node randomly
					randNum = rand.nextInt(i);
					g[i][randNum] = 1;
					g[randNum][i] = 1;					
				}
				
				else {                                           // select destination node proportionaly to its degree
					degCumulative[0] = countElements(g[0]);
					for (k = 1; k < i; k++)
						degCumulative[k] = degCumulative[k-1] + 
								countElements(g[k]);
					len = k-1;
				    randNum = rand.nextInt(degCumulative[len-1] + 1);	
					k = 0;
				    while (randNum > degCumulative[k])
				    	k++;
				    
				    g[i][k] = 1;
					g[k][i] = 1;
				    
				}								
			}			
		}
		
		return g;
	}
	
	
	private static int countElements (byte [] a) {
		int count = 0;
		for (int i = 0; i < a.length; i++)
			if (a[i] == 1) count++;
		return count;
	}
}