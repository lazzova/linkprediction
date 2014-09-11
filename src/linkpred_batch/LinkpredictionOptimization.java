package linkpred_batch;
import java.util.ArrayList;

import org.apache.commons.math3.optim.PointValuePair;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import edu.stanford.nlp.optimization.DiffFunction;
import edu.stanford.nlp.optimization.QNMinimizer;

/**
 * 
 * The main class where everything is being initialized and where optimization is being done
 *
 */
public class LinkpredictionOptimization {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int n = 1000;                                           // number of nodes per graph   
		int f1 = 2;                                             // number of features for the first graph
		int f2 = 3;                                             // numer of features for the second graph
		
		int s = 0;                                              // the node whose links we learn, in this case 0 for each graph
		double alpha = 0;                                       // we use no damping factor within the multiplex, interlayer jump instead
		double interlayer = 0.2;                                // interlayer jump coeffiecient
		double b = 1;//e-6;                                     // WMW function parameter
		double lambda = 1;                                      // regularization parameter 
		double learningRate = 1;                                // learningRate
		
		double [] param = {1, -1, 0.5, -2, 1};                  // parameters vector 2 for the first graph, 3 for the second
		DoubleMatrix1D parameters = new 
				DenseDoubleMatrix1D(param);	
						
		int topN = 10; 
		
		
		
		Network [] graphs = new Network [2];                    // build the graph
		ArtificialGraphGenerator.initialize(f1);
		graphs[0] = (Network) ArtificialGraphGenerator.generate(n, f1, s);
		ArtificialGraphGenerator.initialize(f2);
		graphs[1] = (Network) ArtificialGraphGenerator.generate(n, f2, s);
		
		MultiplexNetwork multiplex = new MultiplexNetwork(graphs, interlayer);
		ArtificialGraphGenerator.buildDandL(multiplex, topN, parameters, alpha);
		
		
		
		QNMinimizer qn = new QNMinimizer(/*15, true*/);
		qn.shutUp();
		qn.useMinPackSearch();
		qn.terminateOnAverageImprovement(false);
		double convergenceTolerance = 1e-10;
		double [] initialGuess = new double [f1+f2];
		int maxFunctionEvaluations = 200;
		LinkPredictionTrainer lp = new LinkPredictionTrainer(
				new RandomWalkGraph [] {multiplex}, f1+f2, alpha, lambda, b, learningRate);
		DiffFunction dfunction = new OptimizationFunction(lp);
		int restarts = 5;
		double [] optimum = null;
		double [] currentOptimum;
		double optimalValue = Double.MAX_VALUE;
		
		GradientDescent gd = new GradientDescent(lp, 200, 1e-6, 6);
		PointValuePair opt = null;
		
		// OPTIMIZATION START
		long start = System.nanoTime();	
		
		while (restarts-- > 0 && optimalValue > 6) {
			for (int i = 0; i < initialGuess.length; i++)
				initialGuess[i] = Math.random() * 2 - 1;
			
			lp.setB(1);
			lp.setLearningRate(1);
			currentOptimum = qn.minimize(dfunction, convergenceTolerance, initialGuess, maxFunctionEvaluations);
			
			lp.setB(1e-6);
			lp.setLearningRate(0.003);
			try {
				opt = gd.optimize(currentOptimum);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			
			
			if (opt.getSecond() < optimalValue) {
				optimalValue = opt.getSecond();
				optimum = opt.getFirst();
			}
		}
		
		long end = System.nanoTime();
		// OPTIMIZATION END
				
		System.out.println();
		System.out.println("\nResults in in " + (end-start)/60E9 + " minutes.");
		System.out.println("Cost: " + optimalValue);
		
		System.out.println();
		System.out.println("TRUE PARAMETERS:");
		for (int i = 0; i < param.length; i++)
			System.out.print(param[i] + " ");
		System.out.println();
		System.out.println("PREDICTED PARAMETERS");
		for (int i = 0; i < optimum.length; i++)
			System.out.print(optimum[i] + " ");
		
		
		// PREDICTIONS
		System.out.println("\n\nPREDICTIONS:");
		ArrayList<Integer> trueLinks = new ArrayList<Integer>();
		ArrayList<Integer> predictedLinks = new ArrayList<Integer>();
				
		Network [] testGraphs = new Network [2];                    // build the graph
		ArtificialGraphGenerator.initialize(f1);
		testGraphs[0] = (Network) ArtificialGraphGenerator.generate(n, f1, s);
		ArtificialGraphGenerator.initialize(f2);
		testGraphs[1] = (Network) ArtificialGraphGenerator.generate(n, f2, s);
				
		MultiplexNetwork testMultiplex = new MultiplexNetwork(graphs, interlayer);
				
		trueLinks = Ranker.predictLinks(testMultiplex, parameters, alpha, topN);
		predictedLinks = Ranker.predictLinks(testMultiplex, new DenseDoubleMatrix1D(optimum), alpha, topN);
				
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


/**
 * 
 * Optimizable function that QNMinimizer uses
 *
 */
class OptimizationFunction implements DiffFunction {
	public LinkPredictionTrainer lp;	

	
	public OptimizationFunction(LinkPredictionTrainer lp) {
		super();
		this.lp = lp;
	}

	@Override
	public double valueAt(double [] x) {
		return lp.getCost(x);		
	}

	@Override
	public int domainDimension() {
		return lp.getParametersNumber();		
	}

	@Override
	public double[] derivativeAt(double [] x) {
		double [] grad = null;
		try {
			grad = lp.getGradient(x);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		
		return grad;
	}
	
}