package teamwork.linkpred_matrix;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import org.apache.commons.math3.util.Pair;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;


class Graph {
	
	int dim;                                                     // number of nodes 
	int s;                                                       // the node whose links we learn
	ArrayList<FeatureField> list;                                // the graph
	ArrayList<Pair<Integer, Double>> D;                          // the future links set
	ArrayList<Pair<Integer, Double>> L;                          // the future no-link set
	SparseCCDoubleMatrix2D A;                                    // the adjacency matrix
	
	// useful
	double [] rowSums;                                           // sum of the each row of the adjacency matrix
	
	
	/**
	 *  Constructor
	 *  
	 *  @param dim
	 */
	public Graph (int n) {
		this.s = 0;
		this.dim = n;
		this.list = new ArrayList<FeatureField>();
		this.A = new SparseCCDoubleMatrix2D(n, n);
		this.D = new ArrayList<Pair<Integer, Double>>();
		this.L = new ArrayList<Pair<Integer, Double>>();
		this.rowSums = new double [dim];                         // filled in the buildTransitionMatrix method
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
		for (int i = 0; i < list.size(); i++) {
			temp = Math.exp(param.zDotProduct(list.get(i).features));
			A.set(list.get(i).row, list.get(i).column, temp);
			A.set(list.get(i).column, list.get(i).row, temp);
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
		SparseCCDoubleMatrix2D Q = buildTransitionMatrix(s, alpha);
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
	public SparseCCDoubleMatrix2D buildTransitionMatrix (int s, double alpha) {
		
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
		// build the transpose of Q  TODO
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
	public DoubleMatrix1D pagerank (SparseCCDoubleMatrix2D Q) {
		
		int n = Q.rows();
		DoubleMatrix1D p = new DenseDoubleMatrix1D(n);           // current iteration
		DoubleMatrix1D oldP = new DenseDoubleMatrix1D(n);        // previous iteration
		//SparseCCDoubleMatrix2D Qtranspose = Q.getTranspose();  TODO
		
		p.assign(1.0 / n);                                       // pagerank initialization 
		
		do {
		
			oldP.assign(p);
			// Qtranspose.zMult(oldP, p);  TODO
			Q.zMult(oldP, p);
					
			oldP.assign(p, new DoubleDoubleFunction() {
		
				@Override
				public double apply(double arg0, double arg1) {
					return Math.abs(arg0-arg1);
				}
			});
		
		} while (oldP.zSum() > 1E-6);                    // convergence check
		
		return p;
	}		
}
	