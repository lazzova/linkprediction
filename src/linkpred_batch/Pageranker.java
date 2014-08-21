package linkpred_batch;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;

public class Pageranker {
	public Graph graph;
	
	// useful
	public double [] rowSums;                                    // sum of the each row of the adjacency matrix
	public DoubleMatrix1D p;                                     // pagerank
	public DoubleMatrix1D [] dp;                                 // pagerank gradient
	
	public Pageranker (Graph graph) {
		this.graph = graph;
		this.rowSums = new double [graph.dim];                         // filled in the buildTransitionMatrix method
		this.p = new DenseDoubleMatrix1D(graph.dim);
		this.dp = new DoubleMatrix1D [graph.f];
		for (int i = 0; i < graph.f; i++)
			dp[i] = new DenseDoubleMatrix1D(graph.dim);
	}
	
	/**
	* Builds the transition matrix for given adjacency matrix
	*
	* @param alpha: damping factor
	* @return SparseCCDoubleMatrix2D
	*/
	public SparseCCDoubleMatrix2D buildTransitionTranspose (double alpha) {
		
		SparseCCDoubleMatrix2D Q = new SparseCCDoubleMatrix2D(graph.A.rows(), graph.A.columns());
		
		// row sums
		int r, c;
		for (int i = 0; i < graph.dim; rowSums[i++] = 0);
		for (int i = 0; i < graph.list.size(); i++) {
			r = graph.list.get(i).row;
			c = graph.list.get(i).column;
			rowSums[r] += graph.A.get(r, c);
			if (r != c)
				rowSums[c] += graph.A.get(c, r);	
		}
		
		// (1-alpha) * A[i][j] / sumElements(A[i])) + 1(j == s) * alpha
		// build the transpose of Q 
		double value;
		for (int i = 0; i < graph.list.size(); i++) {
			r = graph.list.get(i).row;
			c = graph.list.get(i).column;
			value = graph.A.get(r, c);
			value *= (1 - alpha);
			value /= rowSums[r];
			Q.set(c, r, value);
		
			if (r == c) continue;
		
			value = graph.A.get(c, r);
			value *= (1 - alpha);
			value /= rowSums[c];
			Q.set(r, c, value);
		}
		
		for (int i = 0; i < Q.rows(); i++) {
			value = Q.get(graph.s, i);
			value += alpha;
			Q.set(graph.s, i, value);
		}
		
		return Q;				
	}
	
	
	/**
	* Calculates the pagerank, given a transition matrix,
	* using the power method
	* 
	* @param Qt: transpose of the transition probability matrix
	* @return
	*/
	public DoubleMatrix1D pagerank (SparseCCDoubleMatrix2D Qt) {
		
		int n = Qt.rows();
		DoubleMatrix1D p = new DenseDoubleMatrix1D(n);           // current iteration
		DoubleMatrix1D oldP = new DenseDoubleMatrix1D(n);        // previous iteration
				
		p.assign(1.0 / n);                                       // pagerank initialization 
		
		do {
		
			oldP.assign(p);
			Qt.zMult(oldP, p);
					
			oldP.assign(p, new DoubleDoubleFunction() {
		
				@Override
				public double apply(double arg0, double arg1) {
					return Math.abs(arg0-arg1);
				}
			});
		
		} while (oldP.zSum() > 1E-6);                    // convergence check
		
		return p;
	}
	
	
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
	public double edgeWeightPartialD (int nodeIndex, int row, int column, int featureIndex) {		
		return graph.A.get(row, column) * graph.list.get(nodeIndex).features.get(featureIndex);
	}
	
	
	/**
	 * Returns matrix of partial derivatives of the transition matrix
	 *  with respect to the featureIndex-th parameter for the given graph 
     * 
	 * @param featureIndex: the index of the parameter with respect to which the derivative is being calculated 
	 * @param alpha: the damping factor
	 * @return SparseCCDoubleMatrix2D
	 */
	public SparseCCDoubleMatrix2D transitionDerivativeTranspose (int featureIndex, double alpha) {
		
		SparseCCDoubleMatrix2D dQt = new SparseCCDoubleMatrix2D(graph.dim, graph.dim);
		
		// derivative row sums
		int r, c;
		double [] dRowSums = new double [graph.dim];
		for (int i = 0; i < graph.list.size(); i++) {
			r = graph.list.get(i).row;
			c = graph.list.get(i).column;
			dRowSums[r] += edgeWeightPartialD(i, r, c, featureIndex);
			if (r != c)
				dRowSums[c] +=edgeWeightPartialD(i, c, r, featureIndex);	
		}
		
		double value;
		for (int i = 0; i < graph.list.size(); i++) {
			r = graph.list.get(i).row;
			c = graph.list.get(i).column;
			value = (edgeWeightPartialD(i, r, c, featureIndex) * rowSums[r]) -
					(graph.A.get(r, c) * dRowSums[r]);
			value *= (1 - alpha);
			value /= Math.pow(rowSums[r], 2);
			//dQ.set(r, c, value); TODO  Return directly the transpose
			dQt.set(c, r, value);
			
			if (c == r) continue;
			
			value = (edgeWeightPartialD(i, c, r, featureIndex) * rowSums[c]) -
					(graph.A.get(c, r) * dRowSums[c]);
			value *= (1 - alpha);
			value /= Math.pow(rowSums[c], 2);
			//dQ.set(c, r, value);
			dQt.set(r, c, value);
		}
				
		return dQt;
	}
	
	
	/**
	 * Calculates pagerank and it's gradient, for given graph index
	 *  
	 * @param param: the parameters for building the adjacency matrix
	 * @param alpha: the damping factor
	 */
	public void pageRankAndGradient (DoubleMatrix1D param, double alpha) {
		graph.buildAdjacencyMatrix(param);
		SparseCCDoubleMatrix2D Qt = buildTransitionTranspose(alpha);
		
		double EPSILON = 1e-6;
		DoubleMatrix1D oldP = new DenseDoubleMatrix1D(graph.dim);        // the value of p in the previous iteration
						
		DoubleMatrix1D oldDp = new DenseDoubleMatrix1D(graph.dim);       // the value of dp in the previous iteration
		                                                           // ...starts with all entries 0 
		// PAGERANK GRADIENT
		DoubleMatrix1D tmp = new DenseDoubleMatrix1D(graph.dim);;
		for (int k = 0; k < graph.f; k++) {                              // for every parameter
			oldDp.assign(DoubleFunctions.constant(0));
			p.assign(1.0 / graph.dim); 
			dp[k].assign(DoubleFunctions.constant(0));             
			do {
				oldDp.assign(dp[k]);
				
				//transitionDerivative(k, alpha).getTranspose().zMult(p, tmp); TODO
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
