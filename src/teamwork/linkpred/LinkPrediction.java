package teamwork.linkpred;

import java.util.HashMap;

public class LinkPrediction {
	
	static final byte EXPONENTIAL = 1;                           // edge-strength functions 
	static final byte LOGISTIC = 2;
    	
	private ListGraph graph;
	private double alpha;                                        // damping factor
	private double lambda;                                       // regularization parameter
	private double b;                                            // b parameter for the WMW loss function
	private double [] p;                                         // page rank
	private double [][] dp;                                      // page rank gradient
	private double J;                                            // cost
	private double gradient;                                     // gradient of J
	private double [] wPredicted;                                // predicted parameters
	private HashMap<String, Double> q;                           // transition matrix
	// We call it the transition matrix 
    // by convention but we represent it 
	// with a HashMap
	private int s;                                               // index of the s node
	private byte weightFunction;                                 // edge-streingth function (1 or 2)
	private byte [] d;                                           // nodes s will link to in the future 
	                                                             //     if i-th node is in d, d[i] = 1  
	
	
	public LinkPrediction (ListGraph graph) {
		/** Constructor */
		this.graph = graph;
	    this.q = new HashMap<String,Double> ();
	    this.p = new double [graph.n];
	    this.dp = new double [graph.m][graph.n];
	}
	
	public void setS (int sIndex) {
		/** Select the S node, whose links the
		 *  algorithm learns from
		 */
		this.s = sIndex;
	}
	
	public void setWeightFunction (byte weightFunction) {
		/** Sets the edge-strength function */
		this.weightFunction = weightFunction;
	}

	public void edgeWeigth (double [] parameters) {
		/** Sets the weights for all edges in the graph according
		 *  to the chosen edge strength function
		 */
		for (int i = 0; i < graph.n; i++) {
			for (int j = 0; j < graph.adjList[i].size(); j++) {
				if (this.weightFunction == EXPONENTIAL)
					graph.adjList[i].get(j).weight = exponential (dotProduct(
							graph.adjList[i].get(j).features, parameters));
				else if (this.weightFunction == LOGISTIC)
					graph.adjList[i].get(j).weight = logistic (dotProduct(
							graph.adjList[i].get(j).features, parameters));
			}
		}
	}
	
	public static double dotProduct (double [] v1, double [] v2) {
		/** Calculate dot product between two vectors */
		double dProd = 0;
		for (int i = 0; i < v1.length; i++)
			dProd += (v1[i] * v2[i]);
		return dProd;
	}
	
	private double exponential (double z) {
		/** Calculates the exponential function */
		return Math.exp(z);
	}
	
	private double [] exponentialGradient (Edge e, double [] w) {
		/** Calculate partial derivatives of the exponential function
		 *  parameterized by w, with respect to the feature 
		 *  vector e.features
		 */
		double [] der = new double [w.length];
		double tmp = exponential(dotProduct(e.features, w));
		for (int i = 0; i < w.length; i++) {
			der[i] = tmp * w[i];
		}
		return der;
	}
	
	private double logistic (double z) {
		/** Calculates the logistic function */
		return 1.0 / (1+ Math.exp(-z));
	}
	
	private double [] logisticGradient (Edge e, double [] w) {
		/** Calculate partial derivatives of the logistic function
		 *  parameterized by w, with respect to the feature 
		 *  vector e.features
		 */
		double [] der = new double [w.length];
		double tmp = logistic(dotProduct(e.features, w));
		tmp *= (1 - tmp);
		for (int i = 0; i < w.length; i++) {
			der[i] = tmp * w[i];
		}
		return der;
	}
	
	public double [] edgeWeightGradient (Edge e, double [] w) {
		/** Calculate partial derivatives of the weight function
		 *  parameterized by w, with respect to the feature 
		 *  vector e.features
		 */
		if (weightFunction == EXPONENTIAL)
			return exponentialGradient(e, w);
		if (weightFunction == LOGISTIC)
			return logisticGradient(e, w);
		return null;
	}
	
	public void setAlpha(double alpha) {
		/** Sets the damping factor */
		this.alpha = alpha;
	}
	
	public void buildTransitionMatrix () {
		/** Builds the transition matrix q */
		double tmp = 0;
		for (int i = 0; i < graph.n; i++) {
			for (int j = 0; j < graph.adjList[i].size(); i++) {
				tmp = (1-alpha) * graph.adjList[i].get(j).weight / 
						  graph.sumWeights(i);
				if (j == s)
					tmp += alpha;
				q.put(String.format("%d,%d", i, j), tmp);
			}
		}
	}
	
	public double transitionDerivative (Edge e, int index) {
		/** Calculates the partial derivatives of the (u,v) entry
		 *  of the transition matrix with respect to the index-th parameter
		 */
		double derivative = 0;
				
		// TODO
		
		return derivative;
	}
	
	public void calculatePageRank () {
		/** Calculates the page rank, using pagerank
		 *  with restarts algorithm
		 */
		double EPSILON = 1e-12;
		double diff = Double.MAX_VALUE;
		String key;
		double [] oldP = p.clone();
		
		while (diff > EPSILON) {
			diff = 0;
			for (int idx = 0; idx < p.length; idx++) {
				p[idx] = 0;
				for (int i = 0; i < graph.n; i++) {
					key = String.format("%d,%d", i, idx);
					if (q.containsKey(key)) 
						p[idx] += oldP[i] * q.get(key);												
				}
				diff += Math.abs(p[idx] - oldP[idx]);
			}
			oldP = p.clone();
		}				
	}
	
	public void pageRankGradient () {
		/** Calculates the gradient of the page rank */
		for (int i = 0; i < graph.n; i++) 
			p[i] = 1.0 / graph.n;
		
		double EPSILON = 1e-12;                                  // maximum difference allowed between iterations 
		double diff;                                             // current difference between iterations  
		
		double [] oldP = p.clone();                              // the value of p in the previous iteration
		double [][] oldDp = dp.clone();                          // the value of dp in the previous iteration
		
		double tmp;
		double dtk;                                              // k-th partial derivative of the transition matrix
		String key;
		for (int k = 0; k < graph.m; k++) {                      // for every parameter
			diff = Double.MAX_VALUE;
			while (diff > EPSILON) {
				diff = 0;
				for (int u = 0; u < graph.n; u++) {
					tmp = 0;
					for (int v = 0; v < graph.adjList[u].size(); v++) {
						dtk = transitionDerivative(graph.adjList[u].get(v), k);
						tmp += q.get(String.format("%d,%d", v, u)) * 
							   oldDp[k][v] + oldP[v] * dtk;
					}
					dp[k][u] = tmp;
					diff += Math.abs(dp[k][u] - oldDp[k][u]);
					oldDp[k][u] = dp[k][u];
				}
				
				// calculate next iteration page rank
				for (int idx = 0; idx < p.length; idx++) {
					p[idx] = 0;
					for (int i = 0; i < graph.n; i++) {
						key = String.format("%d,%d", i, idx);
						if (q.containsKey(key)) 
							p[idx] += oldP[i] * q.get(key);												
					}
				}
				oldP = p.clone();				
			}		
		}
	}
		
	
	public void buildD (int size) {
		/** Build D, the set of nodes that s will
		 *  link to in the future, with length 'size'.
		 *  D contains the first 'size' number of nodes
		 *  with highest page rank that s does not link to.
		 */
		
		d = new byte[graph.n];
						
		double max = -1;
		int maxIndex = -1;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < graph.n; j++) {
				if (j != s &&
					  p[j] > max &&
					    d[j] == 0 &&
						  !graph.adjList[s].contains(j)) {
					max = this.p[j];
					maxIndex = j;
				}
			}
			this.d[maxIndex] = 1;
			max = -1;
		}		
	}
	
	public void setInitialParameters (double [] initial) {
		/** Set initial parameter values */
		this.wPredicted = initial;
	}
	
	public void setB (double b) {
		/** Sets the parameter for the WMW loss function */
		this.b = b;
	}
	
	public void setLambda (double lambda) {
		/** Sets the regularization parameter */
		this.lambda = lambda;
	}
	
	public double WMWloss (double x) {
		/** Calculates the Wilcoxon-Mann-Whitney loss function */
		return 1.0 / (1+ Math.exp(-x/b));
	}
	
	public double WMWgradient (double x) {
		/** Calculates the gradient of the 
		 *  Wilcoxon-Mann-Whitney loss function 
		 */
		double tmp = logistic(x/b);
		return tmp * (1-tmp) / b;     		
	}
	
	public double costFunction () {
		/** Calculates the fitting error J */
		double regTerm = dotProduct(wPredicted, wPredicted);     // regularization term
		
		double errorTerm = 0;                                    // error term
		for (int i = 0; i < d.length; i++) {
			if (d[i] == 1)                                       // has link
				for (int j = 0; j < graph.n; j++) 
					if (d[j] == 0)                               // no link
						errorTerm += WMWloss(p[j] - p[i]);		
		}
		
		return regTerm + lambda * errorTerm;
	}
	
	public double [] costFunctionGradient () {
		// TODO
		return null;
	}	
}
