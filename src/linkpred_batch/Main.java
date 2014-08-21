package linkpred_batch;

import org.apache.commons.math3.optim.PointValuePair;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;


public class Main {
	
	/**
	 * Main
	 *  
	 * @param args
	 */
	public static void main(String[] args) {
		
		System.out.println("BATCH");                            //TODO 
		System.out.println("Graph generation start");           //TODO
		
		// TODO the optimizator throws an exception with 50 graphs, no exception for 20 or less
		int g = 1;                                             // number of graphs   50
		int n = 10000;                                          // number of nodes    10000
		int f = 2;                                              // number of features 2
		
		int s = 0;                                               // the node whose links we learn, in this case 0 for each graph
		double alpha = 0.2;                                      // damping factor
		double b = 1e-6;                                         // WMW function parameter
		double lambda = 1;                                       // regularization parameter 
		double [] param = {1, -1};                               // parameters vector
		DoubleMatrix1D parameters = new DenseDoubleMatrix1D(param);	
		int topN = 10;
		
		ArtifitialGraphGenerator.initRandom(f);                  // build the graph
		Graph [] graph = new Network [g];
		for (int i = 0; i < g; i++)
			graph[i] = ArtifitialGraphGenerator.generate(n, f, s, topN, parameters, alpha);
		
		System.out.println("Graph generation end");				 //TODO	
		
		
		long start = System.nanoTime();
				
		LinkpredProblem problem = new LinkpredProblem(graph, f, alpha, lambda, b);
		problem.optimize();
		PointValuePair optimum = problem.getOptimum();
		
		long end = System.nanoTime();
		
		System.out.println("Function minimum: " + optimum.getValue() + "\nParameters: " + 
		        optimum.getPoint()[0] + " " + optimum.getPoint()[1]);
		
						
		
		System.out.println("Results in in " + (end-start)/60E9 + " minutes.");
		
	}
}
	
	