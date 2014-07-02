package teamwork.linkpred_matrix;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ml.distance.ManhattanDistance;

public class LinkPrediction {
	private int g;                                               // number of graphs
	private int n;                                               // number of nodes
	private int f;                                               // number of features
	private byte [][][] connected;                               // g matrices, each (i, j)th element is 1 if i and j are linked, 0 otherwise
	private double [][][][] graphs;                              // g graphs for which, each (i,j)th element is a feature vector of the link between i and j
	private double alpha;                                        // damping factor
	private double lambda;                                       // regularization parameter
	private double b;                                            // b parameter for the WMW loss function
	private double [] p;                                         // page rank
	private double [][] dp;                                      // page rank gradient
	private double [][] A;                                       // adjacency matrix
	private RealMatrix Q;                                        // transition matrix	
	private int [] s;                                            // array od indices of all s nodes (MIGHT NOT BE NEEDED)
	private byte [][] D;                                         // nodes s will link to in the future 
	                                                             //     if i-th node is in d, d[i] = 1
	                                                             //     i-th row is the D set for s[i]
	
	private double J;                                            // cost
	private double [] gradient;                                  // gradient
	
	
	public LinkPrediction(double [][][][] graphs, byte [][][] connected, 
			int [] s, byte [][] D, double alpha, double lambda, double b) {
		/** Constructor */
		this.g = graphs.length;
		this.n = graphs[0].length;
		this.f = graphs[0][0][0].length;
		this.graphs = graphs;
		this.connected = connected;
		this.alpha = alpha;
		this.lambda = lambda;
		this.s = s;
		this.b = b;                                              // needed olnly if WMW loss function is used
		this.A = new double [n][n];
		this.Q = new BlockRealMatrix(n, n);
		this.p = new double [n];
		this.dp = new double [f][n];                             // partial deriv of each pagerankvalue are column vectors
		this.D = D;
		this.gradient = new double [f];
		this.J = Double.MAX_VALUE;
	}

	
	private void buildAdjacencyMatrix (int k, RealVector param) {
		/** Builds the adjacency matrix for the k-th graph, using exponential 
		 * edge-strength function, logistic function is the other option.
		 */
		for (int i = 0; i < n; i++) 
			for (int j = 0; j < n; j++) 
				A[i][j] = Math.exp(
						param.dotProduct(new ArrayRealVector(graphs[k][i][j])));
	}
	
	
	private void buildTransitionMatrix (int k) {
		/** Builds the transition matrix for given adjacency 
		 *  matrix and k-th graph 
		 */		
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A.length; j++) {
				if (j == s[k])
					Q.setEntry(i, j, alpha);
				if (connected[k][i][j] == 1) 
					Q.setEntry(i, j, Q.getEntry(i, j) + 
						(1-alpha) * A[i][j] / sumElements(A[i]));
			}
		}		
	}
	
	private double sumElements (double [] a) {
		/** Sums the elements of an array */
		double sum = 0;
		for (int i = 0; i < a.length; i++)
			sum += a[i];
		return sum;
	}
	
	
	public double edgeWeightPartialD (int graph, int i, int j, int index) {
		/** Calculate partial derivative of the weight function (exponential funcion 
		 *  considered) parameterized by w, with respect to the index-th parameter
		 *  for the given graph
		 */
		return A[i][j] * graphs [graph][i][j][index];
	}
	
	
	public RealMatrix transitionDerivative (int graph, int column, int index) {
		/** Returns array of partial derivatives of the (i, column) entries
		 *  of the transition matrix with respect to the index-th parameter
		 *  for the given graph
		 */
		// TODO: Try to find more vectorized approach
		// TODO: Testing
		
		RealMatrix d = new BlockRealMatrix(n, 1); 
		double tmp;
		double sum;
		double sumSquared;
		for (int i = 0; i < n; i++) {
			tmp = edgeWeightPartialD(graph, i, column, index);		
			sum = 0;
			for (int j = 0; j < n; j++)
				sum += A[i][j];
			tmp *= sum;
			
			sumSquared = sum * sum;
			
			sum = 0;
			for (int j = 0; j < n; j++)
				sum += edgeWeightPartialD(graph, i, j, index);
			sum *= A[i][column];
			tmp -= sum;
			
			tmp /= sumSquared;
			tmp *= (1-alpha);
			d.setEntry(i, 0, tmp);			
		}
			
		return d;
	}
	
	
	private void pageRankAndGradient (int graph) {
		/** Calculates pagerank and it's gradient, for given graph */
		// TODO : This method can be optimized
		// TODO
		
		double EPSILON = 1e-12;
		double diff = Double.MAX_VALUE;
		ManhattanDistance manhattan = new ManhattanDistance();
				
		double [] oldP = new double [n];                         // the value of p in the previous iteration
		double [][] oldDp = new double [f][n];                   // the value of dp in the previous iteration
		                                                         // ...starts with all entries 0 
		
		for (int i = 0; i < n; i++)                              // pagerank initialization 
			p[i] = 1.0 / (double)n;
								
		// PAGERANK GRADIENT
		for (int k = 0; k < f; k++) {                            // for every parameter
			diff = Double.MAX_VALUE;
			while (diff > EPSILON) {
				for (int u = 0; u < n; u++) {
					dp[k][u] = Q.getColumnMatrix(u).preMultiply(oldDp[k])[0] +
					      transitionDerivative(graph, u, k).preMultiply(p)[0];
				}
				diff = manhattan.compute(dp[k], oldDp[k]);
				
				for (int u = 0; u < n; u++)
					oldDp[k][u] = dp[k][u];
				
				// calculate next iteration page rank
				p = Q.preMultiply(p);								
			}		
		}
		
		
		// PAGERANK
		while (manhattan.compute(p, oldP) > EPSILON) {
			p = Q.preMultiply(p);
			for (int i = 0; i < n; i++)
				oldP[i] = p[i];
		}
	}
	
	
	public double WMWloss (double x) {
		/** Calculates the Wilcoxon-Mann-Whitney loss function */
		// TODO: Testing
		return 1.0 / (1+ Math.exp(-x/b));
	}
	
	public double WMWderivative (double x) {
		/** Calculates the derivative of the 
		 *  Wilcoxon-Mann-Whitney loss function 
		 */
		// TODO: Testing
		double tmp = 1.0 / (1+ Math.exp(x/b));
		return tmp * (1-tmp) / b;     		
	}
	
	
	public void costFunctionAndGradient (RealVector w) {
		/** Calculates the fitting error J 
		 *  given initial parameter vector 
		 */		
		// TODO: Testing (especially with the gradient calculation)
		
		double regTerm = w.dotProduct(w);                        // regularization term
		double errorTerm = 0;                                    // error term
		
		for (int k = 0; k < g; k++) {                            // for each graph
			buildAdjacencyMatrix(k, w);
			buildTransitionMatrix(k);
			pageRankAndGradient(k);
			
			for (int i = 0; i < D[k].length; i++) {
				if (D[k][i] == 1) {                                                  // has link
					for (int j = 0; j < n; j++) { 
						if (D[k][j] == 0) {                                          // no link
							errorTerm += WMWloss(p[j] - p[i]);
							for (int idx = 0; idx < f; idx++) {                      // for each element of the gradient vector
								gradient[idx] += (WMWderivative(p[j] - p[i]) *       // derivative of the error term
			    			    	    (dp[idx][j] - dp[idx][i]));								
							}
						}
					}
				}
			}
		}
		
	    J = regTerm + lambda * errorTerm;
	    
	    for (int idx = 0; idx < f; idx++) {
	    	gradient[idx] *= lambda;
			gradient[idx] += 2 * w.getEntry(idx);                // derivative of the regularization term
			//gradient[idx] *= 0.01;                               // add learning rate TODO
	    }
	}
	
	
	public double getCost (double []  w) {
		/** Calculates cost function and gradient
		 *  and returns cost function value
		 */
		costFunctionAndGradient(MatrixUtils.createRealVector(w));
		return J;
	}
	
	public double [] getGradient () {
		/** Returns the gradient of the cost function */
		return gradient;
	}
	
	public int getParametersNumber () {
		return f;
	}
}
