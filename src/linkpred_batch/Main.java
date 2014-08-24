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
		int g = 2;                                              // number of graphs   50
		int n = 1000;                                           // number of nodes    10000
		int f = 2;                                              // number of features 2
		
		int s = 0;                                               // the node whose links we learn, in this case 0 for each graph
		double alpha = 0.2;                                      // damping factor
		double b = 1e-6;                                         // WMW function parameter
		double lambda = 1;                                       // regularization parameter 
		double [] param = {1, -1};                               // parameters vector
		DoubleMatrix1D parameters = new DenseDoubleMatrix1D(param);	
		int topN = 10;
		
		ArtifitialGraphGenerator.initialize(f);                  // build the graph
		RandomWalkGraph [] graphs = new Network [g];
		for (int i = 0; i < g; i++)
			graphs[i] = ArtifitialGraphGenerator.generate(n, f, s, topN, parameters, alpha);
		
		System.out.println("Graph generation end");				 //TODO	
		
		long start = System.nanoTime();
		/*		
		LinkpredProblem problem = new LinkpredProblem(graphs, f, alpha, lambda, b);
		problem.optimize();
		PointValuePair optimum = problem.getOptimum();
		*/
		int maxIterations = 500;
		double gradientTreshold = 1e-3;
		double costThreshold = 4;
		double [] initialParameters = new double [f];
		for (int i = 0; i < f; i++)
			initialParameters[i] = Math.random();
		
		GradientDescent gd = new GradientDescent(new LinkPredictionTrainer(graphs, f, alpha, lambda, b), 
				maxIterations, 
				gradientTreshold, 
				costThreshold);
		PointValuePair optimum = null;
		try {
			optimum = gd.optimize(initialParameters);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		long end = System.nanoTime();
		
		System.out.println(gd.getStopReason());
		System.out.println("Function minimum: " + optimum.getValue() + "\nParameters: " + 
		        optimum.getPoint()[0] + " " + optimum.getPoint()[1]);
		
						
		
		System.out.println("Results in in " + (end-start)/60E9 + " minutes.");
		
	}
}
	
	