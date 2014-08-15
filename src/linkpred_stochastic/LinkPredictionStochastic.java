package linkpred_stochastic;

import linkpred_batch.Graph;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;


public class LinkPredictionStochastic {
	int g;                                               // number of graphs
	int f;                                               // number of features
	Graph [] graphs;                                     // all the graphs
	double alpha;                                        // damping factor
	double lambda;                                       // regularization parameter
	double b;                                            // b parameter for the WMW loss function
	double J;                                            // cost
	double [] gradient;                                  // gradient
	DoubleMatrix1D parameters;                           // parameters
	
	
	/**
	 * Constructor
	 * 
	 * @param graphs
	 * @param f
	 * @param s
	 * @param D
	 * @param alpha
	 * @param lambda
	 * @param b
	 */
	public LinkPredictionStochastic(Graph [] graphs, int f, double alpha, double lambda, double b) {
		
		this.g = graphs.length;
		this.f = f;
		this.graphs = graphs;
		this.alpha = alpha;
		this.lambda = lambda;
		this.b = b;                                              // needed olnly if WMW loss function is used
		this.gradient = new double [f];
		this.J = Double.MAX_VALUE;
	}
	
	
	/**
	 * Calculates the Wilcoxon-Mann-Whitney loss function
	 * 
	 * @param x
	 * @return double
	 */
	public double WMWloss (double x) {
		double res = 1.0 / (1+ Math.exp(-x/b));
		return res;
	}
	
	
	/**
	 * Calculates the derivative of the Wilcoxon-Mann-Whitney 
	 * loss function with respect to (pl - pd)
     * 
	 * @param x
	 * @return double
	 */
	public double WMWderivative (double x) {
		double tmp = 1.0 / (1+ Math.exp(x/b));
		double res = tmp * (1-tmp) / b; 
		return res;     		
	}
	
	
	/**
	 * Calculates the fitting error J 
	 * given initial parameter vector 
	 * 
	 * @param w - parameters vector
	 * @param k - graph index
	 * @throws InterruptedException 
	 */
	public void costFunctionAndGradient (DoubleMatrix1D w, int k) throws InterruptedException {			
		this.parameters = w;
		double regTerm = w.zDotProduct(w);                       // regularization term
		double errorTerm = 0;                                    // error term
		
		for (int i = 0; i < f; i++)                              // clear the gradient
			gradient[i] = 0;
		
		Graph graph = graphs[k];
		graph.pageRankAndGradient(parameters, alpha);       // for each graph 
		
		int l, d;
		double delta;                                                                // pl - pd
		for (int i = 0; i < graph.D.size(); i++) {                           // has link
			for (int j = 0; j < graph.L.size(); j++) {                       // no link
				l = graph.L.get(j).getKey();
				d = graph.D.get(i).getKey();
				delta = graph.p.get(l) - graph.p.get(d);
									
				errorTerm += WMWloss(delta);
				
				for (int idx = 0; idx < f; idx++) {                              // for each element of the gradient vector
					gradient[idx] += (WMWderivative(delta) *                     // derivative of the error term
	    		    	    (graph.dp[idx].get(l) - graph.dp[idx].get(d)));					 
				}
			}
		}			
		
		
	    J = regTerm + lambda * errorTerm;
	    
	    for (int idx = 0; idx < f; idx++) {
	    	gradient[idx] *= lambda;
			gradient[idx] += (2 * w.get(idx));                     // derivative of the regularization term				
	    }
	}
	
	
	/**
	 * Calculates cost function and gradient
	 * and returns cost function value
	 * 
	 * @param w
	 * @return double
	 */
	public double getCost (double []  w) {
		return J;
	}
	
	
	/**
	 * Returns the gradient of the cost function
	 * 
	 * @return double []
	 * @throws InterruptedException 
	 */
	public double [] getGradient (double [] w, int k) throws InterruptedException {
		costFunctionAndGradient(new DenseDoubleMatrix1D(w), k);
		return gradient;
	}
	
	
	public int getParametersNumber () {
		return f;
	}
}

