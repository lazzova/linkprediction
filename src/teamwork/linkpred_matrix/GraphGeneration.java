package teamwork.linkpred_matrix;

import java.util.Date;

import org.apache.commons.math3.random.GaussianRandomGenerator;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomVectorGenerator;
import org.apache.commons.math3.random.UncorrelatedRandomVectorGenerator;

public class GraphGeneration {
	
	public static JDKRandomGenerator rand = new JDKRandomGenerator();
	public static RandomVectorGenerator randomVector = null;
	
	/**
	 * Set seed and initialize the generator of
	 * random vectors of length f
	 *
	 * @param f
	 */
	public static void initRandom (int f) {
		rand.setSeed(new Date().getTime());
		randomVector = 
				new UncorrelatedRandomVectorGenerator(
				f, new GaussianRandomGenerator(rand));
	}
	
	/**
	 * Generate graph in FeatureMatrix model
	 * Used for testing purpose	
	 *  
	 * @param n
	 * @return FeatureMatrix
	 */
	public static Graph generate (int n) {
		Graph fm = new Graph(n);
				
		int [] degCumulative = new int [n];	                     // array for cumulative degree sums
		int [] deg = new int [n];                                // array of node degrees
			
		int k;
		int len;
		int randNum;
		                     
		fm.add(0, 1, randomVector.nextVector());          // connect first three nodes in a triad
		fm.add(1, 2, randomVector.nextVector());
		fm.add(2, 0, randomVector.nextVector());
		deg[0] = 2;
		deg[1] = 2;
		deg[2] = 2;
					
		for (int i = 3; i < n; i++) {
			for (int j = 0; j < 3; j++) {                        // generate three links
				if (rand.nextInt(11) < 8) {                      // select destination node randomly
					randNum = rand.nextInt(i);
					fm.add(i, randNum, randomVector.nextVector());
					deg[i]++;
					deg[randNum]++;
				}
				
				else {                                           // select destination node proportionally to its degree
					degCumulative[0] = deg[0];
					for (k = 1; k < i; k++)
						degCumulative[k] = degCumulative[k-1] + deg[k];
					len = k-1;
				    randNum = rand.nextInt(degCumulative[len-1] + 1);	
					k = 0;
				    while (randNum > degCumulative[k])
				    	k++;
				    
				    fm.add(i, k, randomVector.nextVector());
				    deg[i]++;
				    deg[k]++;					    
				}								
			}			
		}		
		
		return fm;
	}
}