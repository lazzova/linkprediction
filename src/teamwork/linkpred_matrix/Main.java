package teamwork.linkpred_matrix;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ml.distance.ManhattanDistance;
import org.apache.commons.math3.optim.PointValuePair;

public class Main {

	
	
	/*
	public static void main(String[] args) {
		// EXPECTED RESULTS:      
		// X1 = 1.5
		// x2 = -1.75
		
		NonLinearConjugateGradientOptimizer opt = new NonLinearConjugateGradientOptimizer(
				NonLinearConjugateGradientOptimizer.Formula.POLAK_RIBIERE, 
				new SimplePointChecker(1.0e-5, 1.0e-3, 100));   
		// create optimizer using Polak-Ribiere formula, and convergence after 10^-5 difference 
		// between iterations, or 10^-3 cost function, or maximum 100 iterations
		
		ObjectiveFunction func = new ObjectiveFunction(new Function());  // Selecting our function as objective
		ObjectiveFunctionGradient grad = new ObjectiveFunctionGradient(new Gradient()); // Selecting self-defined gradient as objective gradient
		
		JDKRandomGenerator rand = new JDKRandomGenerator();
		rand.setSeed(new Date().getTime()); // creates random generator
		RandomVectorGenerator rvg = new UncorrelatedRandomVectorGenerator(
				2, new GaussianRandomGenerator(rand));  // generates random vector of length 2
		
		MultiStartMultivariateOptimizer optimizer = 
				new MultiStartMultivariateOptimizer(opt, 10, rvg);  // creates multistart optimizer
		
		PointValuePair optimum =
	        optimizer.optimize(func, GoalType.MINIMIZE, grad, new InitialGuess(rvg.nextVector()), new MaxEval(100)); 
		// runs the optimization
				        
	    System.out.println("Function minimum: " + optimum.getValue() + "\nParameters: " + 
	        optimum.getPoint()[0] + " " + optimum.getPoint()[1]);
		
		
		/*
		LinkpredOptimizer lp = new LinkpredOptimizer(null);   // the optimization problem
		
		JDKRandomGenerator rand = new JDKRandomGenerator();
		rand.setSeed(new Date().getTime());
		RandomVectorGenerator rvg = new UncorrelatedRandomVectorGenerator(
				2, new GaussianRandomGenerator(rand));
		MultiStartMultivariateOptimizer optimizer = 
				new MultiStartMultivariateOptimizer(lp, 10, rvg);
		
		PointValuePair optimum = optimizer.optimize(optData)
		/
	}
	*/
	
	public static void main(String[] args) {
		
		System.out.println("START");
		
		double [][][][] graphs = new double [1][][][];
		graphs[0] = Graph.generate(100, 2);
		
		System.out.println("Graph generated...");
		
		double alpha = 0.2;
		double b = 1e-2;
		double lambda = 1;
		double [] parameters = {1, -1}; 
		byte [][] D = new byte [1][];
		D[0] = buildD(graphs[0], 0, MatrixUtils.createRealVector(parameters), 10);
		
		System.out.println("D set found...");
		
		LinkpredProblem problem = new LinkpredProblem(graphs, D, alpha, lambda, b);
		
		System.out.println("Optimization started...");
		
		problem.optimize();
		PointValuePair optimum = problem.getOptimum();
		
		System.out.println("Function minimum: " + optimum.getValue() + "\nParameters: " + 
		        optimum.getPoint()[0] + " " + optimum.getPoint()[1]);
	}
	
	
	public static byte [] buildD (double [][][] graph, int s, RealVector trueParameters, int topN) {
		/** Builds the D set (created links) for synthetic graph and known
		 *  parameter values, by taking the first topN highest ranked nodes.
		 *  s is the node whose links we are looking at
		 */
		int n = graph.length;
		byte [] D = new byte [n];
		int [] nodes = new int [n];
		for (int i = 1; i < n; i++)
			nodes[i] = i;
		
		// find pageranks
		double [][] A = buildAdjacencyMatrix(graph, trueParameters);
		RealMatrix Q = buildTransitionMatrix(A);
		double [] rank = pagerank(Q);
		
		// sort the ranks
		double keyRank;
		int keyNode;
		int k;
		for (int i = 1; i < n; i++) {
		    keyRank = rank[i];
		    keyNode = nodes[i];
		    k = i - 1;
		    while (k >= 0 && rank[k] > keyRank) {
		    	rank[k + 1] = rank[k];
		    	nodes[k + 1] = nodes[k];
		        k--;
		    }
		    rank[k + 1] = keyRank;
		    nodes[k + 1] = keyNode;
		}
		
		// find top ranked nodes
		int count = topN;
		k = 0;
		while (k < n && count > 0) {
			if (nodes[k] != s && A[s][nodes[k]] == 1) {          // the node is not s, and has no links to s previously 
				D[nodes[k]] = 1;
				count--;
			}
			k++;
		}
		
		return D;		
	}
	
	private static double [][] buildAdjacencyMatrix (double [][][] graph, RealVector param) {
		/** Builds the adjacency matrix using exponential 
		 * edge-strength function, logistic function is the other option.
		 */
		int n = graph.length;
		double [][] A = new double [n][n];
		for (int i = 0; i < n; i++) 
			for (int j = 0; j < n; j++)
				A[i][j] = Math.exp(
						param.dotProduct(new ArrayRealVector(graph[i][j])));
		
		return A;
	}
	
	
	private static RealMatrix buildTransitionMatrix (double [][] A) {
		/** Builds the transition matrix for given adjacency matrix*/
		RealMatrix Q = new BlockRealMatrix(A.length, A.length);
		
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A.length; j++) {
				Q.setEntry(i, j, A[i][j] / sumElements(A[i]));
			}
		}
		
		return Q;
	}
	
	private static double sumElements (double [] a) {
		/** Sums the elements of an array */
		int sum = 0;
		for (int i = 0; i < a.length; i++)
			sum += a[i];
		return sum;
	}
	
	private static double [] pagerank (RealMatrix Q) {
		/** Calculates the pagerank, given a transition matrix */
		
		double EPSILON = 1e-12;
		ManhattanDistance manhattan = new ManhattanDistance();
		int n = Q.getRowDimension();
		double [] p = new double [n];
		double [] oldP = new double [n];
		
		for (int i = 0; i < n; i++)                      // pagerank initialization 
			p[i] = 1.0 / n;
		
		while (manhattan.compute(p, oldP) > EPSILON) {
			p = Q.preMultiply(p);
			oldP = p.clone();
		}
		
		return p;
	}
}




