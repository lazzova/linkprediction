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
	private static LinkPrediction problem;
	private static Random r;
	
	public static void main(String[] args) {
		/** Main method for defining and solving the learning
		 *  problem as well as validation
		 */
		graphGeneration(100, 2);                                 // generate graph  
		problem = new LinkPrediction(g);                         // create new link prediction problem
		problem.setS(0);                                         // select the s node
		problem.setWeightFunction(LinkPrediction.EXPONENTIAL);   // set the edge-strength function
		double [] w_real = {0.5, -0.3};                          // choose real parameter values
		problem.edgeWeigth(w_real);                              // assign weights (adjacency matrix)
		problem.setAlpha(0.2);                                   // set the damping factor
		problem.buildTransitionMatrix();                         // build the transition matrix
		//problem.calculatePageRank();                             // calculate the page rank 
		//problem.buildD(20);                                      // build d, the set of nodes s will link to in the future   
		
		//problem.setInitialParameters(randomUniformArray(2));     // inital random values for the parameters 
		//problem.setB(1e-6);                                      // set the parameter for the WMW loss function
		//problem.setLambda(1);                                    // set the regularization parameter
		
	}	
	
	
	private static void graphGeneration (int n, int m) {
		/** Generate undirected graph with n nodes and m
		 *  gaussian features and weightFunction as weight function.
		 *  Used for testing purpose
		 *  Adjacency list representation
		 */
		g = new ListGraph(n, m);
		r = new Random(new Date().getTime());
		int [] degCumulative = new int [n];	                     // array for cumulative degree sums	
				
		g.addEdge(0, 1, randomGausianArray(m));                  // connect first three nodes in a triad
		g.addEdge(1, 2, randomGausianArray(m));
		g.addEdge(2, 0, randomGausianArray(m));
		
		int k;
		int len;
		int randNum;
		Edge e;
		for (int i = 3; i < n; i++) {
			
			for (int j = 0; j < 3; j++) {                        // generate three links
				if (r.nextInt(11) < 8) {                         // select destination node randomly
					e = new Edge(i, r.nextInt(i), randomGausianArray(m));										
				}
				
				else {                                           // select destination node proportionaly to its degree
					degCumulative[0] = g.adjList[0].size();
					for (k = 1; k < i; k++)
						degCumulative[k] = degCumulative[k-1] + 
							               g.adjList[k].size();
					len = k-1;
				    randNum = r.nextInt(degCumulative[len-1] + 1);	
					k = 0;
				    while (randNum > degCumulative[k])
				    	k++;
				    
				    e = new Edge(i, k, randomGausianArray(m));				   
				}
				
				if (!g.hasEdge(e)) 
					g.addEdge(e);				
			}
			
		}
	}
	
	
	private static double [] randomGausianArray (int n) {
		/** Generate array of random gaussian distributed 
		 *  numbers, of length n
		 */
		double [] rga = new double [n];
		for (int i = 0; i < n; i++) 
		    rga[i] = r.nextGaussian();	
		
		return rga;		
	}
	
	
	private static double [] randomUniformArray (int n) {
		/** Generate array of random uniform distributed
		 *  numbers, of length n
		 */
		// TODO: Testing
		Random r = new Random(new Date().getTime());
		double [] rua = new double [n];
		for (int i = 0; i < n; i++) 
		    rua[i] = r.nextDouble();	
		
		return rua;
	}	
}
