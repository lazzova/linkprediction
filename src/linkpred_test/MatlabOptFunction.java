package linkpred_test;

import linkpred_batch.ArtificialGraphGenerator;
import linkpred_batch.LinkPredictionTrainer;
import linkpred_batch.MultiplexNetwork;
import linkpred_batch.Network;
import linkpred_batch.RandomWalkGraph;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;

/**
 * 
 * This class will be used with optimization in Matlab
 *
 */
public class MatlabOptFunction {
	private LinkPredictionTrainer problem;
	private double cost;
	private double [] gradient;
	
	
	/**
	 * Constructor
	 */
	public MatlabOptFunction() {}
	
	
	/**
	 * Builds the graphs, sets the parameters and initializes
	 * the LinkPrediction problem.
	 * 
	 * @param g - number of graphs
	 * @param n - number of nodes
	 * @param f - number of features
	 * @param s - index of the source node
	 * @param alpha - damping factor
	 * @param b - WMW loss function's parameter
	 * @param lambda - regularization parameter
	 * @param param - true parameters
	 */
	public void initProblem (int g, int n, int [] fs, int s, 
			double alpha, double b, double lambda, double [] param, 
			double learningRate, double interlayer, int topN) {
		int f = 0;
		for (int i = 0; i < g; i++)
			f += fs[i];
				
		DoubleMatrix1D parameters = new 
				DenseDoubleMatrix1D(param);	
		DoubleMatrix1D [] params = new DoubleMatrix1D [g];
		
		int start = 0;		
		for (int i = 0; i < g; start += fs[i], i++)
			params[i] = parameters.viewPart(start, fs[i]);
			
		Network [] graphs = new Network [g];                    // build the graph
		for (int i = 0; i < g; i++) {
			ArtificialGraphGenerator.initialize(fs[i]);
			graphs[i] = (Network) ArtificialGraphGenerator.generate(n, fs[i], s);
		}
		
				
		MultiplexNetwork multiplex = new MultiplexNetwork(graphs, interlayer);
		ArtificialGraphGenerator.buildDandL(multiplex, topN, parameters, alpha);
				
		System.out.println("Graph generation end");			   
			
		problem = new LinkPredictionTrainer(new RandomWalkGraph [] {multiplex}, f, alpha, lambda, b, learningRate);		
	}

	
	/**
	 * Calculates the cost and the gradient
	 * 
	 * @param w
	 */
	public void runProblem (double [] w) {
		try {
			gradient = problem.getGradient(w);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		cost = problem.getCost(w);		
	}
	
	
	/**
	 * Returns the cost
	 * 
	 * @return double
	 */
	public double getCost() {
		return cost;
	}

	
	/**
	 * Returns the gradeint
	 * 
	 * @return double []
	 */
	public double[] getGradient() {
		return gradient;
	}	
	
}
