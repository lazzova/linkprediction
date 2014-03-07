/** 
 * 
 * @author Verica Lazova, Magdalena Skrizevska, Vesna Dimovska, Blagoja Kojcev
 * @description This is the main class where all the learning and predicting is done
 * 
 * TODO:
 * adjacency_mat
 * compute_pagerank
 * cost_function
 * edge_strength_derivative
 * edge_strength_function
 * estimate_parameters
 * exponential_edge_strength
 * exponential_edge_strength_derivative
 * learn
 * learn_parallel
 * pagerank_simple
 * predict
 * transition_derivative
 * transition_mat
 * WMW_derivatime
 * WMW_loss
 * 
 */

package teamwork.linkpred;

import java.util.Date;
import java.util.Random;


public class Main {
	private static ListGraph g;
	
	// edge-weight function
	public static byte EXPONENTIAL = 1;
	public static byte LOGISTIC = 2;

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double [] w_real = {0.5, -0.3};
		graphGeneration(100, 2, EXPONENTIAL);
	}	
	
	
	private static void graphGeneration (int n, int m, byte weightFunction) {
		/** Generate undirected graph with n nodes and 
		 *  m features and weightFunction as weight function.
		 *  Used for testing purpose
		 *  Adjacency list representation
		 */
		Random r = new Random(new Date().getTime());
		// array for cumulative degree sums
	    int [] degCumulative = new int [n];		
		ListGraph g = new ListGraph(n);
		
		// connect first three nodes in a triad
		g.addEdge(0, 1, weightFunction, randomGausianArray(m));
		g.addEdge(1, 2, weightFunction, randomGausianArray(m));
		g.addEdge(2, 0, weightFunction, randomGausianArray(m));
		
		int k;
		int len;
		int randNum;
		for (int i = 3; i < n; i++) {
			// generate three links
			for (int j = 0; j < 3; j++) {
				// select destination node randomly
				if (r.nextInt() % 11 < 8) {
					g.addEdge(i, r.nextInt() % i, weightFunction, 
						      randomGausianArray(m));
				}
				// select destination node proportionaly to its degree
				else {
					degCumulative[0] = g.getAdjList()[0].size();
					for (k = 1; k < i; k++)
						degCumulative[k] = degCumulative[k-1] + 
							               g.getAdjList()[k].size();
					len = k-1;
				    randNum = r.nextInt() % degCumulative[len-1];	
					k = 0;
				    while (randNum <= degCumulative[k])
				    	k++;
				    g.addEdge(i, k, weightFunction, 
						      randomGausianArray(m));
				}
			}
		}
		
		Main.g = g;
	}
	
	
	private static double [] randomGausianArray (int n) {
		/** Generate array of random gaussian distributed 
		 *  numbers, of length n
		 */
		Random rGausian = new Random(new Date().getTime());
		double [] rga = new double [n];
		for (int i = 0; i < n; i++) 
		    rga[i] = rGausian.nextGaussian();	
		
		return rga;		
	}
	
	
	private static double [] randomUniformArray (int n) {
		/** Generate array of random uniform distributed
		 *  numbers, of length n
		 */
		Random r = new Random(new Date().getTime());
		double [] rua = new double [n];
		for (int i = 0; i < n; i++) 
		    rua[i] = r.nextDouble();	
		
		return rua;
	}
		
	public static double [] pageRank () {
		return null;
	}
	
	
	
	
}
