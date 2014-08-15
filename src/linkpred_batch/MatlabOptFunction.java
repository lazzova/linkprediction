package linkpred_batch;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;

/**
 * 
 * This class will be used with optimization in Matlab
 *
 */
public class MatlabOptFunction {
	private LinkPrediction problem;
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
		GraphGeneration.initRandom(f);                          // build the graph
		Graph [] graph = new Graph [g];
		for (int i = 0; i < g; i++)
			graph[i] = GraphGeneration.generate(n, f);
		
		DoubleMatrix1D parameters = new DenseDoubleMatrix1D(param);
		
		int topRanked = 10;                                      // building D set (linked set)
		for (int i = 0; i < g; i++)
			graph[i].buildD(topRanked, parameters, s, alpha);
		
		problem = new LinkPrediction(graph, f, alpha, lambda, b);		
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
