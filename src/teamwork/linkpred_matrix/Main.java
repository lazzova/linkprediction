package teamwork.linkpred_matrix;

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
		
		System.out.println("Graph generation start");            //TODO
		
		int g = 20;                                             // number of graphs   50
		int n = 100;                                            // number of nodes    10000
		int f = 2;                                              // number of features 2
		
		GraphGeneration.initRandom(f);                          // build the graph
		Graph [] graph = new Graph [g];
		for (int i = 0; i < g; i++)
			graph[i] = GraphGeneration.generate(n, f);
		
		System.out.println("Graph generation end");				 //TODO
		
		int s = 0;                                               // the node whose links we learn, in this case 0 for each graph
		double alpha = 0.2;                                      // damping factor
		double b = 1e-6;                                         // WMW function parameter
		double lambda = 1;                                       // regularization parameter 
		double [] param = {1, -1};                               // parameters vector
		DoubleMatrix1D parameters = new DenseDoubleMatrix1D(param);		
		
		System.out.println("Building D start");                  //TODO
		
		int topRanked = 10;                                      // building D set (linked set)
		for (int i = 0; i < g; i++)
			graph[i].buildD(topRanked, parameters, s, alpha);
				
		long start = System.nanoTime();
		
		
		LinkpredProblem problem = new LinkpredProblem(graph, f, alpha, lambda, b);
		problem.optimize();
		PointValuePair optimum = problem.getOptimum();
		System.out.println("Function minimum: " + optimum.getValue() + "\nParameters: " + 
		        optimum.getPoint()[0] + " " + optimum.getPoint()[1]);
		       
				
		long end = System.nanoTime();
		System.out.println("Results in in " + (end-start)/60E9 + " minutes.");
		
	}
}
	
	