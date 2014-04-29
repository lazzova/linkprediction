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
	
	
	public static double [][][] generate (int n, int f) {
		/** Generate undirected graph with n nodes and m
		 *  Gaussian features
		 *  Used for testing purpose		 
		 */
				
		double [][][] g = new double [n][n][f];		
		
		JDKRandomGenerator rand = new JDKRandomGenerator();
		rand.setSeed(new Date().getTime()); 
		RandomVectorGenerator randomVector = 
				new UncorrelatedRandomVectorGenerator(
				f, new GaussianRandomGenerator(rand));           // random vector generator (f elements)
		
		int [] degCumulative = new int [n];	                     // array for cumulative degree sums
		double rv [];                                            // random vector, temporary storage
				
		rv = randomVector.nextVector();
		g[0][1] = rv.clone();                                    // connect first three nodes in a triad
		g[1][0] = rv.clone();
		rv = randomVector.nextVector();
		g[1][2] = rv.clone();  
		g[2][1] = rv.clone(); 
		rv = randomVector.nextVector();
		g[2][0] = rv.clone(); 
		g[0][2] = rv.clone(); 
		
		int k;
		int len;
		int randNum;
		for (int i = 3; i < n; i++) {
			
			for (int j = 0; j < 3; j++) {                        // generate three links
				if (rand.nextInt(11) < 8) {           // select destination node randomly
					randNum = rand.nextInt(i);
					rv = randomVector.nextVector();
					
					if (g[i][randNum] == null) {
						g[i][randNum] = rv.clone();
						g[randNum][i] = rv.clone();
					}
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
				    
				    rv = randomVector.nextVector();
				    if (g[i][k] == null) {
					    g[i][k] = rv.clone();
					    g[k][i] = rv.clone();
				    }
				}								
			}			
		}
		
		return g;
	}
	
	
	private static int countElements (double [][] a) {
		int count = 0;
		for (int i = 0; i < a.length; i++)
			if (a[i] != null) count++;
		return count;
	}
}