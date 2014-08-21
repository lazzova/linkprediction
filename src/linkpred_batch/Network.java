package linkpred_batch;

import java.util.ArrayList;

import cern.colt.matrix.tdouble.DoubleMatrix1D;

public class Network extends Graph {
	
	/**
	 * Constructor
	 * 
	 * @param n: number of nodes
	 * @param f: number of features
	 * @param s: index of starting node
	 * @param list: the graph given as an array of features
	 * @param D: the linked set
	 * @param L: the no-link set
	 */
	public Network (int n, int f, int s, ArrayList<FeatureField> list, 
			ArrayList<Integer> D, ArrayList<Integer> L) {
		super(n, s, f, list, D, L);		
	}

	
	/**
	 * Constructor
	 */
	public Network () {
		super ();
	}
	
		
	/**
	 * Get the Adjacency matrix of the graph using exponential function
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
			temp = Math.exp(param.zDotProduct(list.get(i).features));
			A.set(r, c, temp);
			if (r != c)
				A.set(c, r, temp);
		}		
	}	
}
	