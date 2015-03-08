package linkpred_test;

import graph_generators.DorogovtsevMendesGraphGenerator;
import graph_generators.GraphGenerator;
import graph_generators.PreferentialAttachmentGraphGenerator;
import graph_generators.RandomEuclideanGraphGenerator;
import graph_generators.RandomGraphGenerator;
import graph_generators.ScaleFreeGraphGenerator;
import graph_generators.SmallWorldGraphGenerator;

import java.util.ArrayList;

import linkpred_batch.ArtificialGraphGenerator;
import linkpred_batch.GradientDescent;
import linkpred_batch.LinkPredictionTrainer;
import linkpred_batch.Network;
import linkpred_batch.RandomWalkGraph;
import linkpred_batch.Ranker;

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
		
		int g = 50;                                                  // number of graphs   
		int n = 1000;                                                // number of nodes per graph
		int f = 2;                                                   // number of features per node
		
		int s = 0;                                                   // the starting node
		double alpha = 0.2;                                          // damping factor
		double b = 1e-6;                                             // WMW function parameter
		double lambda = 1;                                           // regularization parameter 
		double [] param = {1, -1};                                   // parameters vector
		DoubleMatrix1D parameters = new DenseDoubleMatrix1D(param);	
		int topN = 10;
		
		/*
		ArtificialGraphGenerator.initialize(f);                       // build the artificial graph
		RandomWalkGraph [] graphs = new Network [g];
		for (int i = 0; i < g; i++) {
			graphs[i] = ArtificialGraphGenerator.generate(n, f, s);
			ArtificialGraphGenerator.buildDandL(graphs[i], topN, parameters, alpha);  
		}
		*/
		GraphGenerator gg = new SmallWorldGraphGenerator(10, 0.3);
		RandomWalkGraph [] graphs = new Network [g];
		for (int i = 0; i < g; i++) {
			graphs[i] = gg.generate(n, f, s);
			gg.buildDandL(graphs[i], topN, parameters, alpha);  
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
		int restarts = 5;
		double gradientTreshold = 1e-3;                              // Gradient convergence threshold  
		double costThreshold = 15;                                   // Minimal cost
		double [] initialParameters = new double [f];
		for (int i = 0; i < f; i++)
			initialParameters[i] = Math.random();
		
		double learningRate = 0.0001;
		GradientDescent gd = new GradientDescent(
				new LinkPredictionTrainer(graphs, f, alpha, lambda, b, learningRate ), 
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
		//RandomWalkGraph testGraph = ArtificialGraphGenerator.generate(n, f, s);
		RandomWalkGraph testGraph = gg.generate(n, f, s);
		trueLinks = Ranker.predictLinks(
					testGraph, new DenseDoubleMatrix1D(trueParameters), alpha, topN);
		predictedLinks = Ranker.predictLinks(testGraph, new DenseDoubleMatrix1D(optimum.getFirst()), alpha, topN);
				
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
	
	