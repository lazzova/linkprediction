package teamwork.linkpred_matrix;

import java.util.concurrent.ArrayBlockingQueue;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ml.distance.ManhattanDistance;

import teamwork.linkpred.Edge;

public class LinkPrediction {
	private int g;                                               // number of graphs
	private int n;                                               // number of nodes
	private int f;                                               // number of features
	private double [][][][] graphs;                              // g graphs for which, each (i,j)th element is a feature vector of the link between i and j
	private double alpha;                                        // damping factor
	private double lambda;                                       // regularization parameter
	private double b;                                            // b parameter for the WMW loss function
	private double [] p;                                         // page rank
	private double [][] dp;                                      // page rank gradient
	private double [][] A;                                       // adjacency matrix
	private RealMatrix Q;                                        // transition matrix	
	private int [] s;                                            // array od indices of all s nodes
	private byte [][] D;                                         // nodes s will link to in the future 
	                                                             //     if i-th node is in d, d[i] = 1
	                                                             //     i-th row is the D set for s[i]
	
	
	public LinkPrediction(double [][][][] graphs, int [] s, byte [][] D, double alpha, double lambda, double b) {
		/** Constructor */
		this.g = graphs.length;
		this.n = graphs[0].length;
		this.f = graphs[0][0][0].length;
		this.graphs = graphs;
		this.alpha = alpha;
		this.lambda = lambda;
		this.b = b;                                              // needed olnly if WMW loss function is used
		this.s = s;
		this.A = new double [n][n];
		this.Q = new BlockRealMatrix(n, n);
		this.p = new double [n];
		this.dp = new double [f][n];                            // partial deriv of each pagerankvalue are column vectors
		this.D = D;
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
	
	
	private void buildTransitionMatrix () {
		/** Builds the transition matrix for given adjacency matrix*/		
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A.length; j++) {
				Q.setEntry(i, j, A[i][j] / sumElements(A[i]));
			}
		}		
	}
	
	private double sumElements (double [] a) {
		/** Sums the elements of an array */
		int sum = 0;
		for (int i = 0; i < a.length; i++)
			sum += a[i];
		return sum;
	}
	
	
	public double transitionDerivative (int from, int to, int index) {
		/** Calculates the partial derivatives of the (from, to) entry
		 *  of the transition matrix with respect to the index-th parameter
		 */
		// TODO
		
		/*
		Edge e = graph.getEdge(from, to);
		if (e == null)
			return 0;
		
		double derivative = edgeWeightPartialD(e, index);		
		double tmp = 0;
		for (int i = 0; i < graph.adjList[from].size(); i++)
			tmp += graph.adjList[from].get(i).weight;
		derivative *= tmp;
		
		double tmpSquare = tmp * tmp;
		tmp = 0;
		for (int i = 0; i < graph.adjList[from].size(); i++)
			tmp += edgeWeightPartialD(
					graph.adjList[from].get(i), index);
		tmp *= e.weight;
		derivative -= tmp;
		
		derivative /= tmpSquare;
		derivative *= alpha;
		*/	
			
		return 0;
	}
	
	
	private void pageRankAndGradient () {
		/** Calculates pagerank and it's gradient */
		// TODO : This method can be optimized
		// TODO
		
		double EPSILON = 1e-12;
		double diff = Double.MAX_VALUE;
		ManhattanDistance manhattan = new ManhattanDistance();
				
		for (int i = 0; i < n; i++)                              // pagerank initialization 
			p[i] = 1.0 / n;
		
		
		double [] oldP = p.clone();                              // the value of p in the previous iteration
		double [][] oldDp = dp.clone();                          // the value of dp in the previous iteration
		
		double tmp;
		double dtk;                                              // k-th partial derivative of the transition matrix
		String key;
		for (int k = 0; k < f; k++) {                            // for every parameter
			diff = Double.MAX_VALUE;
			while (diff > EPSILON) {
				diff = 0;
				for (int u = 0; u < n; u++) {
					tmp = 0;
					tmp += Q.getColumnMatrix(u).preMultiply(oldDp[k])[0];
					for (int v = 0; v < n; v++) {
						// TODO
						// bi bilo dobro da vrakja matrica (kolona vektor) 
						// nx1 kade sekoj i-ti element kje bide k-ti parcijalen 
						// izvod za Q(i,u),  vo toj slucaj i za ovoj del kje
						// mozime da primenime mnozenje na matrici i da go 
						// izbegneme for ciklusot
						dtk = transitionDerivative(u, v, k);
						tmp += Q.getEntry(v, u) * oldDp[k][v] 
								+ oldP[v] * dtk;
					}
					dp[k][u] = tmp;
					diff += Math.abs(dp[k][u] - oldDp[k][u]);
					oldDp[k][u] = dp[k][u];
				}
				
				// calculate next iteration page rank
				p = Q.preMultiply(p);
				oldP = p.clone();				
			}		
		}
		
		
		// PAGERANK
		while (manhattan.compute(p, oldP) > EPSILON) {
			p = Q.preMultiply(p);
			oldP = p.clone();
		}
	}
	
}
