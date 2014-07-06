package teamwork.linkpred_matrix;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;

public class LinkPrediction {
	private int g;                                               // number of graphs
	private int n;                                               // number of nodes
	private int f;                                               // number of features
	private Graph [] graphs;                                     // all the graphs
	private double alpha;                                        // damping factor
	private double lambda;                                       // regularization parameter
	private double b;                                            // b parameter for the WMW loss function
	private DoubleMatrix1D p;                                    // page rank
	private DoubleMatrix1D [] dp;                                // page rank gradient 
	private SparseCCDoubleMatrix2D Q;                            // transition matrix 
	private int [] s;                                            // array of indices of all s nodes 
	private double J;                                            // cost
	private double [] gradient;                                  // gradient
	
	
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
	public LinkPrediction(Graph [] graphs, int f, int [] s, double alpha, double lambda, double b) {
		
		this.g = graphs.length;
		this.n = graphs[0].dim;
		this.f = f;
		this.graphs = graphs;
		this.alpha = alpha;
		this.lambda = lambda;
		this.s = s;
		this.b = b;                                              // needed olnly if WMW loss function is used
		this.Q = null;
		this.p = new DenseDoubleMatrix1D(n);
		this.dp = new DoubleMatrix1D [f];                        // partial derivative of each pagerank value are column vectors 
		for (int i = 0; i < f; i++)
			dp[i] = new DenseDoubleMatrix1D(n);
		this.gradient = new double [f];
		this.J = Double.MAX_VALUE;
	}

	
	/**
	 * Calculate partial derivative of the weight function (exponential funcion 
	 * considered) parameterized by w, with respect to the index-th parameter
	 * for the given graph
	 * 
	 * @param graphIndex
	 * @param nodeIndex
	 * @param featureIndex
	 * @return double
	 */
	public double edgeWeightPartialD (int graphIndex, int nodeIndex, int row, int column, int featureIndex) {
		
		return graphs[graphIndex].A.get(row, column) * 
			   graphs[graphIndex].list.get(nodeIndex).features.get(featureIndex);
	}
	
	
	/**
	 * Returns matrix of partial derivatives of the transition matrix
	 *  with respect to the featureIndex-th parameter for the given graph 
     * 
	 * @param graph
	 * @param featureIndex
	 * @return SparseCCDoubleMatrix2D
	 */
	public SparseCCDoubleMatrix2D transitionDerivative (int graph, int featureIndex) {
				
		// TODO: Testing
		SparseCCDoubleMatrix2D dQ = new SparseCCDoubleMatrix2D(n, n);
		
		// derivative row sums
		int r, c;
		double [] dRowSums = new double [n];
		for (int i = 0; i < graphs[graph].list.size(); i++) {
			r = graphs[graph].list.get(i).row;
			c = graphs[graph].list.get(i).column;
			dRowSums[r] += edgeWeightPartialD(graph, i, r, c, featureIndex);
			if (r != c)
				dRowSums[c] +=edgeWeightPartialD(graph, i, r, c, featureIndex);	
		}
		
		double value;
		for (int i = 0; i < graphs[graph].list.size(); i++) {
			r = graphs[graph].list.get(i).row;
			c = graphs[graph].list.get(i).column;
			value = (edgeWeightPartialD(graph, i, r, c, featureIndex) * graphs[graph].rowSums[r]) -
					(graphs[graph].A.get(r, c) * dRowSums[r]);
			value *= (1 - alpha);
			value /= Math.pow(graphs[graph].rowSums[r], 2);
			dQ.set(r, c, value);
			
			if (c == r) continue;
			
			value = (edgeWeightPartialD(graph, i, c, r, featureIndex) * graphs[graph].rowSums[c]) -
					(graphs[graph].A.get(c, r) * dRowSums[c]);
			value *= (1 - alpha);
			value /= Math.pow(graphs[graph].rowSums[c], 2);
			dQ.set(c, r, value);
		}
				
		return dQ;
	}
	
	
	/**
	 * Calculates pagerank and it's gradient, for given graph index
	 *  
	 * @param graph
	 */
	private void pageRankAndGradient (int graph) {
				
		// TODO : Testing
		double EPSILON = 1e-6;
		DoubleMatrix1D oldP = new DenseDoubleMatrix1D(n);        // the value of p in the previous iteration
		SparseCCDoubleMatrix2D Qtranspose = Q.getTranspose();  
				
		DoubleMatrix1D oldDp = new DenseDoubleMatrix1D(n);       // the value of dp in the previous iteration
		                                                         // ...starts with all entries 0 
		p.assign(1.0 / n);                                       // pagerank initialization 
				
		// PAGERANK GRADIENT
		DoubleMatrix1D tmp = new DenseDoubleMatrix1D(n);;
		for (int k = 0; k < f; k++) {                            // for every parameter
			oldDp.assign(DoubleFunctions.constant(0));
			p.assign(1.0 / n);                                   // TODO not sure if new initialization is needed each time
			do {
				oldDp.assign(dp[k]);
				
				transitionDerivative(graph, k).getTranspose().zMult(p, tmp); 
				Qtranspose.zMult(oldDp, dp[k]); 
				dp[k].assign(tmp, DoubleFunctions.plus);
				
				oldDp.assign(dp[k], new DoubleDoubleFunction() {
					
					@Override
					public double apply(double arg0, double arg1) {
						return Math.abs(arg0-arg1);
					}
				});
				
				// calculate next iteration page rank
				Qtranspose.zMult(p.copy(), p);  	
				
			} while (oldDp.zSum() > EPSILON);		
		}
		
		// PAGERANK
		do {
			
			//for (int i = 0; i < n; i++)                        TODO
			//	oldP.set(i, p.get(i));
			oldP.assign(p);
			Qtranspose.zMult(oldP, p); 
								
			oldP.assign(p, new DoubleDoubleFunction() {
		
				@Override
				public double apply(double arg0, double arg1) {
					return Math.abs(arg0-arg1);
				}
			});
		
		} while (oldP.zSum() > EPSILON);                         // convergence check
	}
	
	
	/**
	 * Calculates the Wilcoxon-Mann-Whitney loss function
	 * 
	 * @param x
	 * @return double
	 */
	public double WMWloss (double x) {
		return 1.0 / (1+ Math.exp(-x/b));
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
		return tmp * (1-tmp) / b;     		
	}
	
	
	/**
	 * Calculates the fitting error J 
	 * given initial parameter vector 
	 * 
	 * @param w
	 */
	public void costFunctionAndGradient (DoubleMatrix1D w) {			
		// TODO: Testing (especially with the gradient calculation)
		
		double regTerm = w.zDotProduct(w);                       // regularization term
		double errorTerm = 0;                                    // error term
		
		for (int k = 0; k < g; k++) {                            // for each graph
			graphs[k].buildAdjacencyMatrix(w);
			Q = graphs[k].buildTransitionMatrix(s[k], alpha);
			pageRankAndGradient(k); 
		
     		int l, d;
			double delta;                                                            // pl - pd
			for (int i = 0; i < graphs[k].D.size(); i++) {                           // has link
				for (int j = 0; j < graphs[k].L.size(); j++) {                       // no link
					l = graphs[k].L.get(j).getKey();
					d = graphs[k].D.get(i).getKey();
					delta = p.get(l) - p.get(d);
										
					errorTerm += WMWloss(delta);
					for (int idx = 0; idx < f; idx++) {                              // for each element of the gradient vector
						gradient[idx] += (WMWderivative(delta) *                     // derivative of the error term
	    			    	    (dp[idx].get(l) - dp[idx].get(d)));					 
					}
				}
			}			
		}
		
	    J = regTerm + lambda * errorTerm;
	    
	    for (int idx = 0; idx < f; idx++) {
	    	gradient[idx] *= lambda;
			gradient[idx] += 2 * w.get(idx);                       // derivative of the regularization term
			//gradient[idx] *= 0.01;                               // add learning rate TODO
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
		costFunctionAndGradient(new DenseDoubleMatrix1D(w));
		return J;
	}
	
	
	/**
	 * Returns the gradient of the cost function
	 * 
	 * @return double []
	 */
	public double [] getGradient () {
		return gradient;
	}
	
	
	public int getParametersNumber () {
		return f;
	}
}
