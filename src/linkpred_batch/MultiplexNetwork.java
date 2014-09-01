package linkpred_batch;

import java.util.ArrayList;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;

public class MultiplexNetwork extends RandomWalkGraph {
	//public int dim;                                              // number of nodes 
	//public int s;                                                // the node whose links we learn
	//public int f;                                                // number of features
	//public ArrayList<FeatureField> list;                         // the graph
	//public ArrayList<Integer> D;                                 // the future links set
	//public ArrayList<Integer> L;                                 // the future no-link set
	//public SparseCCDoubleMatrix2D A;                             // the adjacency matrix

	public int graphsNumber;                                       // the number of graphs within the multiplex
	public int layerDim;                                           // dimension of a single layer
	public Network [] graphs;
	public double interlayer;                                      // interlayer transition coefficient
	
	// useful
	public double [] rowSums; 
	
	public MultiplexNetwork (Network [] graphs, double interlayer) {
		super();
		this.graphsNumber = graphs.length;
		this.graphs = graphs;
		this.interlayer = interlayer;		
		this.layerDim = graphs[0].dim;
				
		this.dim = graphsNumber * layerDim;
		this.rowSums = new double [dim];
		this.s = graphs[0].s;                                      // s node for the first layer TODO
		for (int i = 0; i < graphsNumber; i++) 
			this.f += graphs[i].f;
		
		this.list = new ArrayList<FeatureField>();
		list.addAll(graphs[0].list);
		for (int i = 1; i < graphsNumber; i++) 
			for (FeatureField f : graphs[i].list) 
				list.add(new FeatureField(i*layerDim+f.row, i*layerDim+f.column, f.features));		
				
		this.D = new ArrayList<Integer>();
		this.L = new ArrayList<Integer>();
		for (int i = 1; i < graphsNumber; i++) {
			for (Integer nodeIndex : graphs[i].D) 
				this.D.add(i * layerDim + nodeIndex);
			for (Integer nodeIndex : graphs[i].L) 
				this.L.add(i * layerDim + nodeIndex);			
		}
		
		this.A = new SparseCCDoubleMatrix2D(dim, dim);
		this.p = new DenseDoubleMatrix1D(this.dim);
		this.dp = new DoubleMatrix1D [this.f];
		for (int i = 0; i < this.f; i++)
			dp[i] = new DenseDoubleMatrix1D(this.dim);
	}
		
	/*	
	@Override
	public void buildAdjacencyMatrix(DoubleMatrix1D param) {
		// intralayer
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
		
		// interlayer
		int startI = 0, startJ = 0;
		for (int i = 0; i < this.interlayer.length; i++) {
			for (int j = 0; j < this.interlayer[0].length; j++) {
				if (i != j) {
					startI = i * this.layerDim;
					startJ = j * this.layerDim;
					for (int k = 0; k < this.layerDim; k++)
						A.set(startI+k, startJ+k, this.interlayer[i][j]);
				}				
			}
		}		
	}
	*/

	/**
	 * Build the Adjacency matrix of the graph
	 * 
	 * @param param: the parameters used for building the adjacency matrix
	 * @return
	 */
	@Override
	public void buildAdjacencyMatrix(DoubleMatrix1D param) {
		
		DoubleMatrix1D [] params = new DoubleMatrix1D [graphsNumber];
		for (int i = 0, start = 0; i < graphsNumber; start += graphs[i++].f)
			params[i] = param.viewPart(start, graphs[i].f);
		
		double temp;
		int r, c;
		for (int i = 0; i < list.size(); i++) {
			r = list.get(i).row;
			c = list.get(i).column;
			temp = weightingFunction(params[r / layerDim].zDotProduct(list.get(i).features));
			A.set(r, c, temp);
			if (r != c)
				A.set(c, r, temp);
		}				
	}
	
	
	//TODO: check this method, likley a bug is present
	@Override
	public SparseCCDoubleMatrix2D buildTransitionTranspose(double alpha) {
		SparseCCDoubleMatrix2D Q = new SparseCCDoubleMatrix2D(dim, dim);
		
		// row sums
		int r, c;                                                 // graph, row, column
		for (int i = 0; i < this.dim; rowSums[i++] = 0);
		for (int i = 0; i < this.list.size(); i++) {
			r = this.list.get(i).row;
			c = this.list.get(i).column;
			rowSums[r] += this.A.get(r, c);
			if (r != c)
				rowSums[c] += this.A.get(c, r);	
		}
		
		// (1-alpha-interlayer) * A[i][j] / sumElements(A[i])) + 1(j == s) * alpha
		// build the transpose of Q 
		double value;
		double factors = alpha + interlayer;
		for (int i = 0; i < this.list.size(); i++) {
			r = this.list.get(i).row;
			c = this.list.get(i).column;
			value = this.A.get(r, c);
			value *= (1 - factors);
			value /= rowSums[r];
			Q.set(c, r, value);
		
			if (r == c) continue;
		
			value = this.A.get(c, r);
			value *= (1 - factors);
			value /= rowSums[c];
			Q.set(r, c, value);
		}
		
		// add damping factor 
		for (int i = 0; i < graphsNumber; i++) {
			for (int k = 0; k < layerDim; k++) {
					value = Q.get(i * layerDim + this.s, i * layerDim + k);
					value += alpha;
					Q.set(i * layerDim + this.s, i * layerDim + k, value);				
			}
		}
		
		
		// add interlayer jumps 
		for (int i = 0; i < graphsNumber; i++) {
			for (int j = 0; j < graphsNumber; j++) {
				if (i == j) continue;
				for (int k = 0; k < layerDim; k++) 	
					Q.set(k + j*layerDim, i*layerDim + k, interlayer);				
			}
		}
		
		return Q;		
	}

	
	// TODO: DEBUG
	public void printMatrix (SparseCCDoubleMatrix2D mat) {
		for (int i = 0; i < mat.rows(); i++) {
			for (int j = 0; j < mat.columns(); j++)
				System.out.printf("%5.2f ", mat.get(i, j));
			System.out.println();
		}
		System.out.println();
	}
	
	// TODO: DEBUG
	public void isColumnStochastic (SparseCCDoubleMatrix2D mat) {
		System.out.println();
		for (int i = 0; i < mat.columns(); i++) 
			System.out.printf("%.2f\n", mat.viewColumn(i).zSum());		
	}
	

	/**
	 * Returns matrix of partial derivatives of the transition matrix
	 *  with respect to the featureIndex-th parameter for the given graph 
     * 
	 * @param featureIndex: the index of the parameter with respect to which the derivative is being calculated 
	 * @param alpha: the damping factor
	 * @return SparseCCDoubleMatrix2D
	 */
	@Override
	public SparseCCDoubleMatrix2D transitionDerivativeTranspose(
			int featureIndex, double alpha) {
		// TODO
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
			value *= (1 - alpha - interlayer); // TODO
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
