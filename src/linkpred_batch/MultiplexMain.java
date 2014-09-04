package linkpred_batch;

import java.util.ArrayList;

import org.apache.commons.math3.optim.PointValuePair;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;

public class MultiplexMain {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("Graph generation start");           
		
		int g = 1;                                              // number of graphs   
		int n = 100;                                            // number of nodes per graph   
		int f1 = 2;                                             // number of features for the first graph
		int f2 = 3;                                             // numer of features for the second graph
		
		int s = 0;                                              // the node whose links we learn, in this case 0 for each graph
		double alpha = 0.1;                                     // damping factor
		double interlayer = 0.2;                                // interlayer jump coeffiecient
		double b = 1e-6;                                        // WMW function parameter
		double lambda = 1;                                      // regularization parameter 
		
		double [] param = {1, -1, 0.5, -2, 1};                  // parameters vector 2 for the first graph, 3 for the second
		DoubleMatrix1D parameters = new 
				DenseDoubleMatrix1D(param);	
		DoubleMatrix1D parameters1 = new 
				DenseDoubleMatrix1D(new double [] {1, -1});	
		DoubleMatrix1D parameters2 = new 
				DenseDoubleMatrix1D(new double [] {0.5, -2, 1});	
				
		int topN = 5;
		
		Network [] graphs = new Network [2];                    // build the graph
		ArtificialGraphGenerator.initialize(f1);
		graphs[0] = (Network) ArtificialGraphGenerator.generate(n, f1, s, topN, parameters1, alpha + interlayer);
		ArtificialGraphGenerator.initialize(f2);
		graphs[1] = (Network) ArtificialGraphGenerator.generate(n, f2, s, topN, parameters2, alpha + interlayer);
		
		MultiplexNetwork multiplex = new MultiplexNetwork(graphs, interlayer);
		
				
		System.out.println("Graph generation end");			   
		
		long start = System.nanoTime();		
		
		// GRADIENT DESCENT OPTIMIZATION START
		
		int maxIterations = 300;
		int restarts = 20;
		double gradientTreshold = 1e-5;
		double costThreshold = 3;
		double [] initialParameters = new double [f1 + f2];
		for (int i = 0; i < f1 + f2; i++)
			initialParameters[i] = Math.random();
		
		GradientDescent gd = new GradientDescent(new LinkPredictionTrainer(
				new RandomWalkGraph [] {multiplex}, f1+f2, alpha, lambda, b), 
				maxIterations, 
				gradientTreshold, 
				costThreshold);
		
		PointValuePair optimum = null;
		try {
			//optimum = gd.optimize(initialParameters);
			optimum = gd.multiStartOptimize(restarts, initialParameters);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
				
		// GRADIENT DESCENT OPTIMIZATION END
		
		
		/*
		LinkpredProblem problem = new LinkpredProblem(new RandomWalkGraph [] {multiplex}, f1+f2, alpha, lambda, b);
		problem.optimize();
		PointValuePair optimum = problem.getOptimum();
		*/
		
		long end = System.nanoTime();
		
		System.out.println(gd.getStopReason());
		System.out.println("Function minimum: " + optimum.getValue() + "\nParameters: ");
		for (int i = 0; i < f1+f2; i++)
		        System.out.print(optimum.getPoint()[i] + " ");
		System.out.println("\nResults in in " + (end-start)/60E9 + " minutes.");
		
		
		// PREDICTIONS
		ArrayList<Integer> trueLinks = new ArrayList<Integer>();
		ArrayList<Integer> predictedLinks = new ArrayList<Integer>();
		double [][] trueParameters = new double [2][];
		trueParameters[0] = parameters1.toArray();
		trueParameters[1] = parameters2.toArray();
		for (int i = 0; i < graphs.length; i++)
			trueLinks.addAll(Ranker.predictLinks(
					graphs[i], new DenseDoubleMatrix1D(trueParameters[i]), alpha + interlayer, topN));
		predictedLinks = Ranker.predictLinks(multiplex, new DenseDoubleMatrix1D(optimum.getFirst()), alpha, topN*graphs.length);
		
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