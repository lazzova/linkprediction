package linkpred_test;

import java.util.ArrayList;

import linkpred_batch.ArtificialGraphGenerator;
import linkpred_batch.GradientDescent;
import linkpred_batch.LinkPredictionTrainer;
import linkpred_batch.MultiplexNetwork;
import linkpred_batch.Network;
import linkpred_batch.RandomWalkGraph;
import linkpred_batch.Ranker;

import org.apache.commons.math3.optim.PointValuePair;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;

public class MultiplexMain {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("Graph generation start");           
		
		int n = 1000;                                           // number of nodes per graph   
		int f1 = 2;                                             // number of features for the first graph
		int f2 = 3;                                             // number of features for the second graph
		int f3 = 2;                                             // number of features for the third graph
		
		int s = 0;                                              // the node whose links we learn, in this case 0 for each graph
		double alpha = 0.0;                                     // we use no damping factor within the multiplex, interlayer jump instead
		double interlayer = 0.2;                                // interlayer jump coeffiecient
		double b = 1e-6;                                        // WMW function parameter
		double lambda = 1;                                      // regularization parameter 
		
		double [] param = {1, -1, 0.5, -2, 1, -2, 2};           // parameters vector 2 for the first graph, 3 for the second
		DoubleMatrix1D parameters = new 
				DenseDoubleMatrix1D(param);	
						
		int topN = 10; 
			
		Network [] graphs = new Network [3];                    // build the graph
		ArtificialGraphGenerator.initialize(f1);
		graphs[0] = (Network) ArtificialGraphGenerator.generate(n, f1, s);
		ArtificialGraphGenerator.initialize(f2);
		graphs[1] = (Network) ArtificialGraphGenerator.generate(n, f2, s);
		ArtificialGraphGenerator.initialize(f3);
		graphs[2] = (Network) ArtificialGraphGenerator.generate(n, f3, s); 
		
		MultiplexNetwork multiplex = new MultiplexNetwork(graphs, interlayer);
		ArtificialGraphGenerator.buildDandL(multiplex, topN, parameters, alpha);
		
		// TODO: debugging
		//multiplex.buildAdjacencyMatrix(new DenseDoubleMatrix1D(param));
		//multiplex.printMatrix(multiplex.buildTransitionTranspose(alpha));
		//multiplex.isColumnStochastic(multiplex.buildTransitionTranspose(alpha));
		
		System.out.println("Graph generation end");			   
		
		long start = System.nanoTime();		
		
		// GRADIENT DESCENT OPTIMIZATION START
		
		int maxIterations = 150;
		int restarts = 10;
		double gradientTreshold = 1e-5;
		double costThreshold = 5.6;
		int f = f1+f2+f3;
		double [] initialParameters = new double [f];
		for (int i = 0; i < f; i++)
			initialParameters[i] = Math.random() * 2 - 1;                 // select random value between -1 and 1 instead between 0 and 1   
		
		
		GradientDescent gd = new GradientDescent(new LinkPredictionTrainer(
				new RandomWalkGraph [] {multiplex}, f, alpha, lambda, b, 0.0003), // TODO 
				maxIterations, 
				gradientTreshold, 
				costThreshold);
		PointValuePair optimum = null;
		
		try {
			optimum = gd.multiStartOptimize(restarts);
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
		for (int i = 0; i < f; i++)
		        System.out.print(optimum.getPoint()[i] + " ");
		System.out.println("\nResults in in " + (end-start)/60E9 + " minutes.");
		
				
		// PREDICTIONS
		ArrayList<Integer> trueLinks = new ArrayList<Integer>();
		ArrayList<Integer> predictedLinks = new ArrayList<Integer>();
		
		Network [] testGraphs = new Network [3];                    // build the graph
		ArtificialGraphGenerator.initialize(f1);
		testGraphs[0] = (Network) ArtificialGraphGenerator.generate(n, f1, s);
		ArtificialGraphGenerator.initialize(f2);
		testGraphs[1] = (Network) ArtificialGraphGenerator.generate(n, f2, s);
		ArtificialGraphGenerator.initialize(f3);
		testGraphs[2] = (Network) ArtificialGraphGenerator.generate(n, f3, s);
		
		MultiplexNetwork testMultiplex = new MultiplexNetwork(graphs, interlayer);
		
		trueLinks = Ranker.predictLinks(testMultiplex, parameters, alpha, topN);
		predictedLinks = Ranker.predictLinks(testMultiplex, new DenseDoubleMatrix1D(optimum.getFirst()), alpha, topN);
		
		System.out.println("\nTrue links:");
		for (int i = 0; i < trueLinks.size(); i++)
			System.out.print(trueLinks.get(i) + " ");
		System.out.println();
		
		System.out.println("\nPredicted links:");
		for (int i = 0; i < predictedLinks.size(); i++)
			System.out.print(predictedLinks.get(i) + " ");
		System.out.println();
		
	}
	
}
