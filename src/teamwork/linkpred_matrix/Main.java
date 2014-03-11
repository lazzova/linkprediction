package teamwork.linkpred_matrix;


//UNTESTED

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


import java.util.Date;
import java.util.Random;


public class Main {
	private static MatrixGraph g;
	private static LinkPrediction problem;

	public static void main(String[] args) {
		/** Main method for defining and solving the learning
		 *  problem as well as validation
		 */
		graphGeneration(100, 2);                                 // generate graph  
		problem.setWeightFunction(LinkPrediction.EXPONENTIAL);   // set the edge-strength function
		double [] w_real = {0.5, -0.3};                          // choose real parameter values
		problem.setEdgeWeigths(w_real);                          // assign weights (adjacency matrix)
		problem.setAlpha(0.2);                                   // set the damping factor
		problem.buildTransitionMatrix();                         // build the transition matrix
		problem.calculatePageRank();                             // TODO: calculate the page rank 

	}	


	private static void graphGeneration (int n, int m) {
		/** Generate undirected graph with n nodes and m
		 *  gaussian features and weightFunction as weight function.
		 *  Used for testing purpose
		 *  Adjacency list representation
		 */
		Random r = new Random(new Date().getTime());
		int [] degCumulative = new int [n];	                     // array for cumulative degree sums	

		g.addEdge(0, 1, randomGausianArray(m));                  // connect first three nodes in a triad
		g.addEdge(1, 2, randomGausianArray(m));
		g.addEdge(2, 0, randomGausianArray(m));

		int k;
		int len;
		int randNum;
		for (int i = 3; i < n; i++) {

			for (int j = 0; j < 3; j++) {                        // generate three links
				if (r.nextInt() % 11 < 8) {                      // select destination node randomly
					g.addEdge(i, r.nextInt() % i, randomGausianArray(m));
				}

				else {                                           // select destination node proportionaly to its degree
					degCumulative[0] = g.n;
					for (k = 1; k < i; k++)
						degCumulative[k] = degCumulative[k-1] + 
							               g.adjMatrix[k].length;
					len = k-1;
				    randNum = r.nextInt() % degCumulative[len-1];	
					k = 0;
				    while (randNum <= degCumulative[k])
				    	k++;
				    g.addEdge(i, k, randomGausianArray(m));
				}
			}

		}
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
}
