package linkpred_test;

import java.util.Date;

import linkpred_batch.LinkPredictionTrainer;
import linkpred_batch.RandomWalkGraph;

import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimplePointChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultiStartMultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;
import org.apache.commons.math3.optim.nonlinear.scalar.gradient.NonLinearConjugateGradientOptimizer;
import org.apache.commons.math3.random.GaussianRandomGenerator;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomVectorGenerator;
import org.apache.commons.math3.random.UncorrelatedRandomVectorGenerator;

public class LinkpredProblem {
	/**The trainer object*/
	private LinkPredictionTrainer lp;
	/**The optimal point and cost function value at it*/
	private PointValuePair optimum;
	
	
	/**
	 * Constructor (Creates the trainer object)
	 * 
	 * @param graphs: the graphs
	 * @param f: the number of features per node
	 * @param alpha: damping factor
	 * @param lambda: regularization parameter
	 * @param b: the b parameter of the WMW loss function
	 */
	public LinkpredProblem (RandomWalkGraph [] graphs, int f, double alpha, double lambda, double b) {
		this.lp = new LinkPredictionTrainer(graphs, f, alpha, lambda, b, 0.0003);
	}

	
	/**
	 *  Does the optimization
	 */
	public void optimize () {
		NonLinearConjugateGradientOptimizer opt = new NonLinearConjugateGradientOptimizer(
				NonLinearConjugateGradientOptimizer.Formula.FLETCHER_REEVES, 
				new SimplePointChecker<PointValuePair>(1e-5, 1e-5, 300));                               // create optimizer using FLETCHER-REEVES formula, and convergence 
		                                                                                                // after 10^-4 relative threshold, or 10^4 absolute threshold, 
		                                                                                                // or maximum 100 iterations
		LinkPredictionCost costFunction = new LinkPredictionCost(lp);
		ObjectiveFunction func = new ObjectiveFunction(costFunction);                                   // Selecting self-defined function as objective
		ObjectiveFunctionGradient grad = new ObjectiveFunctionGradient(new LinkPredictionGradient(lp)); // Selecting self-defined gradient as objective gradient
		
		JDKRandomGenerator rand = new JDKRandomGenerator();
		rand.setSeed(new Date().getTime());                                                              // creates random generator
		RandomVectorGenerator rvg = new UncorrelatedRandomVectorGenerator(
				lp.getParametersNumber(), new GaussianRandomGenerator(rand));                            // generates random vector of initial parameters
		
		MultiStartMultivariateOptimizer optimizer = 
				new MultiStartMultivariateOptimizer(opt, 5, rvg);                                        // creates multistart optimizer with 20 starting points
		
		System.out.println("And we are running ...");
		
		optimum = optimizer.optimize(
				func, GoalType.MINIMIZE, grad, new InitialGuess(rvg.nextVector()), new MaxEval(300));    // runs the optimization    
		
		
		System.out.println("\n\n");                                                                      // prints information about the optimization (for debugging purpose)
		System.out.println("Number of evaluations: " + optimizer.getEvaluations());
		System.out.println("Max number of evaluations: " + optimizer.getMaxEvaluations());
		System.out.println("Number of iterations: " + optimizer.getIterations());
		System.out.println("Max number of iterations: " + optimizer.getMaxIterations());
			
	}
	
	
	/**
	 * Returns the optimal point
	 * 
	 * @return PointValuePair
	 */
	public PointValuePair getOptimum() {
		return optimum;
	}
	
}
