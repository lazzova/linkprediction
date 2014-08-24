package linkpred_batch;

import java.util.ArrayList;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
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
		this.rowSums = new double [dim];
		
		this.dim = graphsNumber * layerDim;
		this.s = graphs[0].s;                                      // s node for the first layer TODO
		this.f = graphs[0].f;
		
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
	
	
	@Override
	public SparseCCDoubleMatrix2D buildTransitionTranspose(double alpha) {
		SparseCCDoubleMatrix2D Q = new SparseCCDoubleMatrix2D(this.A.rows(), this.A.columns());
		
		// row sums
		int g, r, c;                                                 // graph, row, column
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
		// TODO
		return Q;		
	}

	
	public void printMatrix (SparseCCDoubleMatrix2D mat) {
		for (int i = 0; i < mat.rows(); i++) {
			for (int j = 0; j < mat.columns(); j++)
				System.out.printf("%5.2f ", mat.get(i, j));
			System.out.println();
		}
		System.out.println();
	}
	

	@Override
	public SparseCCDoubleMatrix2D transitionDerivativeTranspose(
			int featureIndex, double alpha) {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public double weightingFunction(double x) {
		return Math.exp(x);
	}


	@Override
	public double weightingFunctionDerivative(int nodeIndex, int row,
			int column, int featureIndex) {
		// TODO Auto-generated method stub
		return 0;
	}
	

}
