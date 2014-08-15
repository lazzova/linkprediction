package linkpred_batch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import org.apache.commons.math3.util.Pair;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;


public class Graph {
	
	public int dim;                                              // number of nodes 
	public int s;                                                // the node whose links we learn
	public int f;                                                // number of features
	public ArrayList<FeatureField> list;                         // the graph
	public ArrayList<Pair<Integer, Double>> D;                   // the future links set
	public ArrayList<Pair<Integer, Double>> L;                   // the future no-link set
	public SparseCCDoubleMatrix2D A;                             // the adjacency matrix
	
	// useful
	public double [] rowSums;                                    // sum of the each row of the adjacency matrix
    public DoubleMatrix1D p;                                     // pagerank
    public DoubleMatrix1D [] dp;                                 // pagerank gradient
	
	/**
	 *  Constructor
	 *  
	 *  @param dim
	 */
	public Graph (int n, int f) {
		this.dim = n;
		this.f = f;
		this.list = new ArrayList<FeatureField>();
		this.A = new SparseCCDoubleMatrix2D(n, n);
		this.D = new ArrayList<Pair<Integer, Double>>();
		this.L = new ArrayList<Pair<Integer, Double>>();
		this.rowSums = new double [dim];                         // filled in the buildTransitionMatrix method
		this.p = new DenseDoubleMatrix1D(n);
		this.dp = new DoubleMatrix1D [f];
		for (int i = 0; i < f; i++)
			dp[i] = new DenseDoubleMatrix1D(dim);
	}

	
	/**
	 * Add new matrix element
	 * 
	 * @param row
	 * @param column
	 * @param features
	 */
	public void add (int row, int column, double [] features) {
		FeatureField ff = new FeatureField(row, column, features);
		if (!list.contains(ff))
			list.add(new FeatureField(row, column, features));
	}
	
	
	/**
	 * Get the Adjacency matrix of the graph using exponential function
	 * 
	 * @param param
	 * @return
	 */
	public void buildAdjacencyMatrix (DoubleMatrix1D param) {
				
		double temp;
		int r, c;
		for (int i = 0; i < list.size(); i++) {
			r = list.get(i).row;
			c = list.get(i).column;
			temp = Math.exp(param.zDotProduct(list.get(i).features));
			A.set(r, c, temp);
			if (r != c)
				A.set(c, r, temp);
		}		
	}
	
	
	/**
	 * Builds the D set (created links) for synthetic graph and known
	 * parameter values, by taking the first topN highest ranked nodes.
	 * s is the node whose links we are looking at.
	 *
	 * @param topN
	 * @param trueParameters
	 * @param s
	 * @param alpha
	 */
	public void buildD (int topN, DoubleMatrix1D trueParameters, int s, double alpha) {
		this.s = s;
		
		// find pageranks
		buildAdjacencyMatrix(trueParameters);
		SparseCCDoubleMatrix2D Q = buildTransitionTranspose(s, alpha);
		DoubleMatrix1D rank = pagerank(Q);
		
		// sort the ranks in ascending order
		for (int i = 0; i < rank.size(); i++) 
			if (i != s && A.get(s, i) == 0)                      // the node is not s, and has no links to s previously
				L.add(new Pair<Integer, Double> (i, rank.get(i)));
				
		Collections.sort(L, new Comparator<Pair<Integer, Double>> () {
		
			@Override
			public int compare(Pair<Integer, Double> o1,
					Pair<Integer, Double> o2) {
				if (o1.getValue() > o2.getValue()) return 1;
				if (o2.getValue() > o1.getValue()) return -1;
				return 0;
			}		
		});
		
		// put the highest ranked in D and remove those from L
		while (D.size() < topN) {
			D.add(L.get(L.size()-1));
			L.remove(L.size()-1);			
		}		
	}
		
		
	/**
	* Builds the transition matrix for given adjacency matrix
	* and s as starting node
	*
	* @param s
	* @param alpha
	* @return SparseCCDoubleMatrix2D
	*/
	public SparseCCDoubleMatrix2D buildTransitionTranspose (int s, double alpha) {
		
		SparseCCDoubleMatrix2D Q = new SparseCCDoubleMatrix2D(A.rows(), A.columns());
		
		// row sums
		int r, c;
		for (int i = 0; i < dim; rowSums[i++] = 0);
		for (int i = 0; i < list.size(); i++) {
			r = list.get(i).row;
			c = list.get(i).column;
			rowSums[r] += A.get(r, c);
			if (r != c)
				rowSums[c] += A.get(c, r);	
		}
		
		// (1-alpha) * A[i][j] / sumElements(A[i])) + 1(j == s) * alpha
		// build the transpose of Q 
		double value;
		for (int i = 0; i < list.size(); i++) {
			r = list.get(i).row;
			c = list.get(i).column;
			value = A.get(r, c);
			value *= (1 - alpha);
			value /= rowSums[r];
			Q.set(c, r, value);
		
			if (r == c) continue;
		
			value = A.get(c, r);
			value *= (1 - alpha);
			value /= rowSums[c];
			Q.set(r, c, value);
		}
		
		for (int i = 0; i < Q.rows(); i++) {
			value = Q.get(s, i);
			value += alpha;
			Q.set(s, i, value);
		}
		
		return Q;				
	}
	
	
	/**
	* Calculates the pagerank, given a transition matrix,
	* using the power method
	* 
	* @param Q
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
	 * @param graphIndex
	 * @param nodeIndex
	 * @param featureIndex
	 * @return double
	 */
	public double edgeWeightPartialD (int nodeIndex, int row, int column, int featureIndex) {		
		return A.get(row, column) * list.get(nodeIndex).features.get(featureIndex);
	}
	
	
	/**
	 * Returns matrix of partial derivatives of the transition matrix
	 *  with respect to the featureIndex-th parameter for the given graph 
     * 
	 * @param graph
	 * @param featureIndex
	 * @return SparseCCDoubleMatrix2D
	 */
	public SparseCCDoubleMatrix2D transitionDerivativeTranspose (int featureIndex, double alpha) {
		
		SparseCCDoubleMatrix2D dQt = new SparseCCDoubleMatrix2D(dim, dim);
		
		// derivative row sums
		int r, c;
		double [] dRowSums = new double [dim];
		for (int i = 0; i < list.size(); i++) {
			r = list.get(i).row;
			c = list.get(i).column;
			dRowSums[r] += edgeWeightPartialD(i, r, c, featureIndex);
			if (r != c)
				dRowSums[c] +=edgeWeightPartialD(i, c, r, featureIndex);	
		}
		
		double value;
		for (int i = 0; i < list.size(); i++) {
			r = list.get(i).row;
			c = list.get(i).column;
			value = (edgeWeightPartialD(i, r, c, featureIndex) * rowSums[r]) -
					(A.get(r, c) * dRowSums[r]);
			value *= (1 - alpha);
			value /= Math.pow(rowSums[r], 2);
			//dQ.set(r, c, value); TODO  Return directly the transpose
			dQt.set(c, r, value);
			
			if (c == r) continue;
			
			value = (edgeWeightPartialD(i, c, r, featureIndex) * rowSums[c]) -
					(A.get(c, r) * dRowSums[c]);
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
	 * @param graph
	 */
	public void pageRankAndGradient (DoubleMatrix1D param, double alpha) {
		buildAdjacencyMatrix(param);
		SparseCCDoubleMatrix2D Qt = buildTransitionTranspose(s, alpha);
		
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
	