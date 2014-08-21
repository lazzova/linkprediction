package linkpred_batch;

import java.util.ArrayList;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;

public abstract class Graph {
	public int dim;                                              // number of nodes 
	public int s;                                                // the node whose links we learn
	public int f;                                                // number of features
	public ArrayList<FeatureField> list;                         // the graph
	public ArrayList<Integer> D;                                 // the future links set
	public ArrayList<Integer> L;                                 // the future no-link set
	public SparseCCDoubleMatrix2D A;                             // the adjacency matrix
	
	
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
	public Graph(int dim, int s, int f, ArrayList<FeatureField> list,
			ArrayList<Integer> D, ArrayList<Integer> L) {
		this.dim = dim;
		this.s = s;
		this.f = f;
		this.list = list;
		this.D = D;
		this.L = L;
		this.A = new SparseCCDoubleMatrix2D(dim, dim);
	}
	
	
	/**
	 * Constructor
	 */
	public Graph () {
		this.dim = 0;
		this.s = 0;
		this.f = 0;
		this.list = new ArrayList<FeatureField>();
		this.D = new ArrayList<Integer>();
		this.L = new ArrayList<Integer>();
		this.A = null;
	};


	/**
	 * Build the adjacency matrix of the graph given parameters
	 * 
	 * @param param: the parameters used for building the adjacency matrix
	 */
	public abstract void buildAdjacencyMatrix (DoubleMatrix1D param); 
}
