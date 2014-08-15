package teamwork.linkpred_matrix;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;

public class LinkPrediction {
	private int g;                                               // number of graphs
	private int f;                                               // number of features
	private Graph [] graphs;                                     // all the graphs
	private double alpha;                                        // damping factor
	private double lambda;                                       // regularization parameter
	private double b;                                            // b parameter for the WMW loss function
	private double J;                                            // cost
	private double [] gradient;                                  // gradient
	private DoubleMatrix1D parameters;                           // parameters
	
	
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
	public LinkPrediction(Graph [] graphs, int f, double alpha, double lambda, double b) {
		
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
	 * @param w
	 * @throws InterruptedException 
	 */
	public void costFunctionAndGradient (DoubleMatrix1D w) throws InterruptedException {			
		this.parameters = w;
		double regTerm = w.zDotProduct(w);                       // regularization term
		double errorTerm = 0;                                    // error term
		
		for (int i = 0; i < f; i++)                              // clear the gradient
			gradient[i] = 0;
		
		int num_threads = Runtime.getRuntime().availableProcessors()+1;             // TODO concurrency
		ThreadPoolExecutor executor = new ThreadPoolExecutor(num_threads,
				g, Long.MAX_VALUE, TimeUnit.MINUTES, new ArrayBlockingQueue<Runnable>(g));

		for (int k = 0; k < g; k++) { 
			final Graph tmpg = graphs[k];
		    executor.execute(new Runnable() {
									
				@Override
				public void run() {
					try {
						tmpg.pageRankAndGradient(parameters, alpha);       // for each graph 
					}
		            catch (Exception e) {
						e.printStackTrace();
					}
				}			
		    });
		}

		executor.shutdown();
		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.MINUTES);
		
		int l, d;
		double delta;                                            // pl - pd
		for (int k = 0; k < g; k++) {      		                                                            
			for (int i = 0; i < graphs[k].D.size(); i++) {                           // has link
				for (int j = 0; j < graphs[k].L.size(); j++) {                       // no link
					l = graphs[k].L.get(j).getKey();
					d = graphs[k].D.get(i).getKey();
					delta = graphs[k].p.get(l) - graphs[k].p.get(d);
										
					errorTerm += WMWloss(delta);
					
					for (int idx = 0; idx < f; idx++) {                              // for each element of the gradient vector
						gradient[idx] += (WMWderivative(delta) *                     // derivative of the error term
	    			    	    (graphs[k].dp[idx].get(l) - graphs[k].dp[idx].get(d)));					 
					}
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
	public double [] getGradient (double [] w) throws InterruptedException {
		costFunctionAndGradient(new DenseDoubleMatrix1D(w));
		return gradient;
	}
	
	
	public int getParametersNumber () {
		return f;
	}
}
