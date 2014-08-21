package linkpred_batch;

import java.util.ArrayList;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;

public class MultiplexNetwork extends Graph {
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
	public double [][] interlayer;                                 // interlayer links
	
	
	public MultiplexNetwork (Network [] graphs, double [][] interlayer, int graphS, int nodeS,
			ArrayList<Integer> D, ArrayList<Integer> L) {
		super();
		this.graphsNumber = graphs.length;
		this.graphs = graphs;
		this.interlayer = interlayer;
		
		this.layerDim = graphs[0].dim;
		this.dim = graphs.length * graphs[0].dim;
		this.s = graphS * graphs[0].dim;
		this.f = graphs[0].f;
		
		this.list = new ArrayList<FeatureField>();
		list.addAll(graphs[0].list);
		for (int i = 1; i < graphsNumber; i++) {
			for (FeatureField f : graphs[i].list) {
				f.column += i * graphs[i].dim;
				f.row += i * graphs[i].dim;
			}
			this.list.addAll(graphs[i].list);
		}
		
		this.D = D;
		this.L = L;
		this.A = new SparseCCDoubleMatrix2D(dim, dim);
	}
		
		
	@Override
	public void buildAdjacencyMatrix(DoubleMatrix1D param) {
		// intralayer
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
	

}
