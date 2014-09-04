package linkpred_batch;

import java.util.ArrayList;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;

public abstract class RandomWalkGraph {
	/**Number of nodes */
	public int dim;                                              
	/**The starting node*/
	public int s;                                         
	/**Number of features per node*/
	public int f;                                                
	/**The graph*/
	public ArrayList<FeatureField> list;                         
	/**The future links set*/
	public ArrayList<Integer> D;     
	/**The future no-link set*/
	public ArrayList<Integer> L; 
	/**the adjacency matrix*/
	public SparseCCDoubleMatrix2D A;                             
	
	// useful
	/**Pagerank*/
	public DoubleMatrix1D p;                                     
	/**Pagerank gradient*/
	public DoubleMatrix1D [] dp;                                 
	
	
	/**
	 * Constructor
	 * 
	 * @param dim
	 * @param s
	 * @param f
	 * @param list
	 * @param D
	 * @param L
	 */
	public RandomWalkGraph(int dim, int s, int f, ArrayList<FeatureField> list,
			ArrayList<Integer> D, ArrayList<Integer> L) {
		this.dim = dim;
		this.s = s;
		this.f = f;
		this.list = list;
		this.D = D;
		this.L = L;
		this.A = new SparseCCDoubleMatrix2D(dim, dim);
		this.p = new DenseDoubleMatrix1D(this.dim);
		this.dp = new DoubleMatrix1D [this.f];
		for (int i = 0; i < this.f; i++)
			dp[i] = new DenseDoubleMatrix1D(this.dim);
	}
	
	
	/**
	 * Constructor
	 */
	public RandomWalkGraph () {};


	/**
	 * Build the adjacency matrix of the graph given parameters
	 * 
	 * @param param: the parameters used for building the adjacency matrix
	 */
	public abstract void buildAdjacencyMatrix (DoubleMatrix1D param); 
	
	
	/**
	* Build the transition matrix for given adjacency matrix
	*
	* @param alpha: damping factor
	* @return SparseCCDoubleMatrix2D
	*/
	public abstract SparseCCDoubleMatrix2D buildTransitionTranspose (double alpha);
	
	
	/**
	 * Returns matrix of partial derivatives of the transition matrix
	 *  with respect to the featureIndex-th parameter for the given graph 
     * 
	 * @param featureIndex: the index of the parameter with respect to which the derivative is being calculated 
	 * @param alpha: the damping factor
	 * @return SparseCCDoubleMatrix2D
	 */
	public abstract SparseCCDoubleMatrix2D transitionDerivativeTranspose (int featureIndex, double alpha);
	
	
	/**
	 * Defines the edge-weighting function
	 * 
	 * @param x: weighting function argument
	 * @return double
	 */
	public abstract double weightingFunction (double x);
	
	
	/**
	 * Calculate partial derivative of the weight function (exponential funcion 
	 * considered) parameterized by w, with respect to the index-th parameter
	 * for the given graph
	 * 
	 * @param nodeIndex: the index of the node in the graph
	 * @param row: the row index of the adjacency matrix
	 * @param column: the column index of the adjacency matrix
	 * @param featureIndex: the index of the parameter with respect to which the derivative is being calculated 
	 * @return double
	 */
	public abstract double weightingFunctionDerivative (int nodeIndex, int row, int column, int featureIndex);
	
	
	/**
	 * Returns true if a link from 'from' node to 'to' node in the graph,
	 * otherwise returns false
	 * 
	 * @param from: link start node
	 * @param to: link end node
	 * @return boolean
	 */
	public abstract boolean hasLink (int from, int to);
	
	
	/**
	 * Calculates pagerank and it's gradient, for given graph index
	 *  
	 * @param param: the parameters for building the adjacency matrix
	 * @param alpha: the damping factor
	 */
	public void pageRankAndGradient (DoubleMatrix1D param, double alpha) {
		buildAdjacencyMatrix(param);
		SparseCCDoubleMatrix2D Qt = buildTransitionTranspose(alpha);
		
		double EPSILON = 1e-6;
		DoubleMatrix1D oldP = new DenseDoubleMatrix1D(dim);        // the value of p in the previous iteration
						
		DoubleMatrix1D oldDp = new DenseDoubleMatrix1D(dim);       // the value of dp in the previous iteration
		                                                           // ...starts with all entries 0 
		// PAGERANK GRADIENT
		DoubleMatrix1D tmp = new DenseDoubleMatrix1D(dim);;
		for (int k = 0; k < f; k++) {                              // for every parameter
			oldDp.assign(DoubleFunctions.constant(0));
			p.assign(1.0 / dim); 
			dp[k].assign(DoubleFunctions.constant(0));             
			do {
				oldDp.assign(dp[k]);
				
				//transitionDerivative(k, alpha).getTranspose().zMult(p, tmp); 
				transitionDerivativeTranspose(k, alpha).zMult(p, tmp);
				Qt.zMult(oldDp, dp[k]);
				dp[k].assign(tmp, DoubleFunctions.plus);
				
				oldDp.assign(dp[k], new DoubleDoubleFunction() {
					
					@Override
					public double apply(double arg0, double arg1) {
						return Math.abs(arg0-arg1);
					}
				});
				
				// calculate next iteration page rank
				Qt.zMult(p.copy(), p);
			} while (oldDp.zSum() > EPSILON);		
		}
		
		// PAGERANK
		do {
			
			oldP.assign(p);
			Qt.zMult(oldP, p);
								
			oldP.assign(p, new DoubleDoubleFunction() {
		
				@Override
				public double apply(double arg0, double arg1) {
					return Math.abs(arg0-arg1);
				}
			});
		
		} while (oldP.zSum() > EPSILON);                         // convergence check
	}	
	
	
}
