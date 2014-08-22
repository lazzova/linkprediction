package linkpred_batch;

import java.util.ArrayList;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;

public class Network extends RandomWalkGraph {
	// useful
	public double [] rowSums;                                    // sum of the each row of the adjacency matrix
		
	/**
	 * Constructor
	 * 
	 * @param n: number of nodes
	 * @param f: number of features
	 * @param s: index of starting node
	 * @param list: the this given as an array of features
	 * @param D: the linked set
	 * @param L: the no-link set
	 */
	public Network (int n, int f, int s, ArrayList<FeatureField> list, 
			ArrayList<Integer> D, ArrayList<Integer> L) {
		super(n, s, f, list, D, L);		
		this.rowSums = new double [this.dim];                         // filled in the buildTransitionMatrix method
	}

	
	/**
	 * Get the Adjacency matrix of the this using exponential function
	 * 
	 * @param param: the parameters used for building the adjacency matrix
	 * @return
	 */
	public void buildAdjacencyMatrix (DoubleMatrix1D param) {
				
		double temp;
		int r, c;
		for (int i = 0; i < list.size(); i++) {
			r = list.get(i).row;
			c = list.get(i).column;
			temp = weightingFunction(param.zDotProduct(list.get(i).features));
			A.set(r, c, temp);
			if (r != c)
				A.set(c, r, temp);
		}		
	}
	
	
	/**
	* Builds the transition matrix for given adjacency matrix
	*
	* @param alpha: damping factor
	* @return SparseCCDoubleMatrix2D
	*/
	public SparseCCDoubleMatrix2D buildTransitionTranspose (double alpha) {
		
		SparseCCDoubleMatrix2D Q = new SparseCCDoubleMatrix2D(this.A.rows(), this.A.columns());
		
		// row sums
		int r, c;
		for (int i = 0; i < this.dim; rowSums[i++] = 0);
		for (int i = 0; i < this.list.size(); i++) {
			r = this.list.get(i).row;
			c = this.list.get(i).column;
			rowSums[r] += this.A.get(r, c);
			if (r != c)
				rowSums[c] += this.A.get(c, r);	
		}
		
		// (1-alpha) * A[i][j] / sumElements(A[i])) + 1(j == s) * alpha
		// build the transpose of Q 
		double value;
		for (int i = 0; i < this.list.size(); i++) {
			r = this.list.get(i).row;
			c = this.list.get(i).column;
			value = this.A.get(r, c);
			value *= (1 - alpha);
			value /= rowSums[r];
			Q.set(c, r, value);
		
			if (r == c) continue;
		
			value = this.A.get(c, r);
			value *= (1 - alpha);
			value /= rowSums[c];
			Q.set(r, c, value);
		}
		
		for (int i = 0; i < Q.rows(); i++) {
			value = Q.get(this.s, i);
			value += alpha;
			Q.set(this.s, i, value);
		}
		
		return Q;				
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
		
		SparseCCDoubleMatrix2D dQt = new SparseCCDoubleMatrix2D(this.dim, this.dim);
		
		// derivative row sums
		int r, c;
		double [] dRowSums = new double [this.dim];
		for (int i = 0; i < this.list.size(); i++) {
			r = this.list.get(i).row;
			c = this.list.get(i).column;
			dRowSums[r] += weightingFunctionDerivative(i, r, c, featureIndex);
			if (r != c)
				dRowSums[c] +=weightingFunctionDerivative(i, c, r, featureIndex);	
		}
		
		double value;
		for (int i = 0; i < this.list.size(); i++) {
			r = this.list.get(i).row;
			c = this.list.get(i).column;
			value = (weightingFunctionDerivative(i, r, c, featureIndex) * rowSums[r]) -
					(this.A.get(r, c) * dRowSums[r]);
			value *= (1 - alpha);
			value /= Math.pow(rowSums[r], 2);
			//dQ.set(r, c, value); TODO  Return directly the transpose
			dQt.set(c, r, value);
			
			if (c == r) continue;
			
			value = (weightingFunctionDerivative(i, c, r, featureIndex) * rowSums[c]) -
					(this.A.get(c, r) * dRowSums[c]);
			value *= (1 - alpha);
			value /= Math.pow(rowSums[c], 2);
			//dQ.set(c, r, value);
			dQt.set(r, c, value);
		}
				
		return dQt;
	}


	@Override
	public double weightingFunction(double x) {
		return Math.exp(x);
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
	@Override
	public double weightingFunctionDerivative(int nodeIndex, int row, int column, int featureIndex) {
		return this.A.get(row, column) * this.list.get(nodeIndex).features.get(featureIndex);
	}
}
	