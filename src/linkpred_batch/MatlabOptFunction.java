package linkpred_batch;

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
	public void initProblem (int g, int n, int f, int s, 
			double alpha, double b, double lambda, double [] param) {
		ArtificialGraphGenerator.initialize(f);                          // build the graph
		RandomWalkGraph [] graph = new Network [g];
		int topN = 10; 
		DoubleMatrix1D parameters = new DenseDoubleMatrix1D(param);
		
		for (int i = 0; i < g; i++)
			graph[i] = ArtificialGraphGenerator.generate(n, f, s, topN, parameters, alpha);
	
		
		problem = new LinkPredictionTrainer(graph, f, alpha, lambda, b);		
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
