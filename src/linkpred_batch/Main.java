package linkpred_batch;

import java.util.ArrayList;

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
		
		System.out.println("Graph generation start");           
		
		int g = 2;                                                   // number of graphs   
		int n = 100;                                                 // number of nodes per graph
		int f = 2;                                                   // number of features per node
		
		int s = 0;                                                   // the starting node
		double alpha = 0.2;                                          // damping factor
		double b = 1e-6;                                             // WMW function parameter
		double lambda = 1;                                           // regularization parameter 
		double [] param = {0.5, -0.2};                                   // parameters vector
		DoubleMatrix1D parameters = new DenseDoubleMatrix1D(param);	
		int topN = 10;
		
		ArtificialGraphGenerator.initialize(f);                       // build the artificial graph
		RandomWalkGraph [] graphs = new Network [g];
		for (int i = 0; i < g; i++) {
			graphs[i] = ArtificialGraphGenerator.generate(
					n, f, s, parameters, alpha);
			ArtificialGraphGenerator.buildDandL(graphs[i], topN, parameters, alpha);  
		}
		
		System.out.println("Graph generation end");				 
		
		long start = System.nanoTime();
		
		/*		
		LinkpredProblem problem = new LinkpredProblem(graphs, f, alpha, lambda, b);
		problem.optimize();
		PointValuePair optimum = problem.getOptimum();
		*/
		
		// GRADIENT DESCENT OPTIMIZATION START
		
		int maxIterations = 150;                                     // Maximum number of iterations  
		int restarts = 20;
		double gradientTreshold = 1e-3;                              // Gradient convergence threshold  
		double costThreshold = 4;                                    // Minimal cost
		double [] initialParameters = new double [f];
		for (int i = 0; i < f; i++)
			initialParameters[i] = Math.random();
		
		GradientDescent gd = new GradientDescent(
				new LinkPredictionTrainer(graphs, f, alpha, lambda, b), 
				maxIterations, 
				gradientTreshold, 
				costThreshold);
		PointValuePair optimum = null;
		try {
			optimum = gd.multiStartOptimize(restarts);
			//optimum = gd.optimize(initialParameters);                // do the optimization
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		
		// GRADIENT DESCENT OPTIMIZATION END
		
		long end = System.nanoTime();
		
		//System.out.println(gd.getStopReason());
		System.out.println("Function minimum: " + optimum.getValue() + "\nParameters: " + 
		        optimum.getPoint()[0] + " " + optimum.getPoint()[1]);
		
						
		
		System.out.println("Results in in " + (end-start)/60E9 + " minutes.");
		
		
		// PREDICTIONS
		ArrayList<Integer> trueLinks = new ArrayList<Integer>();
		ArrayList<Integer> predictedLinks = new ArrayList<Integer>();
		double [] trueParameters = param;
		RandomWalkGraph testGraph = ArtificialGraphGenerator.generate(n, f, s, parameters, alpha);
		trueLinks = Ranker.predictLinks(
					testGraph, new DenseDoubleMatrix1D(trueParameters), alpha, topN);
		predictedLinks = Ranker.predictLinks(testGraph, new DenseDoubleMatrix1D(optimum.getFirst()), alpha, topN);
				
		System.out.println("\nTrue links:");
		for (int i = 0; i < trueLinks.size(); i++)
		System.out.print(trueLinks.get(i) + " ");
		System.out.println();
		
		System.out.println("\nPredicted links:");
		for (int i = 0; i < predictedLinks.size(); i++)
			System.out.print(predictedLinks.get(i) + "(" + predictedLinks.get(i) % n + ") ");
		System.out.println();
		
	}
}
	
	