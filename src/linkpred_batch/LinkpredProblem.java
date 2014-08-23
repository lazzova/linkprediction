package linkpred_batch;

import java.util.Date;

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
	private LinkPredictionTrainer lp;
	private PointValuePair optimum;
	
	public LinkpredProblem (RandomWalkGraph [] graphs, int f, double alpha, double lambda, double b) {
		this.lp = new LinkPredictionTrainer(graphs, f, alpha, lambda, b);
	}

	public void optimize () {
		NonLinearConjugateGradientOptimizer opt = new NonLinearConjugateGradientOptimizer(
				NonLinearConjugateGradientOptimizer.Formula.FLETCHER_REEVES, 
				new SimplePointChecker<PointValuePair>(0.005, 0.005, 300));                                     // create optimizer using Polak-Ribiere formula, and convergence 
		                                                                                                // after 10^-6 difference between iterations, or 10^2 cost function, 
		                                                                                                // or maximum 100 iterations
		LinkPredictionCost costFunction = new LinkPredictionCost(lp);
		ObjectiveFunction func = new ObjectiveFunction(costFunction);                                   // Selecting our function as objective
		ObjectiveFunctionGradient grad = new ObjectiveFunctionGradient(new LinkPredictionGradient(lp)); // Selecting self-defined gradient as objective gradient
		
		JDKRandomGenerator rand = new JDKRandomGenerator();
		rand.setSeed(new Date().getTime());                                                              // creates random generator
		RandomVectorGenerator rvg = new UncorrelatedRandomVectorGenerator(
				lp.getParametersNumber(), new GaussianRandomGenerator(rand));                            // generates random vector of initial parameters
		
		//MultiStartMultivariateOptimizer optimizer = 
			//	new MultiStartMultivariateOptimizer(opt, 5, rvg);                                        // creates multistart optimizer with 500 starting points
		
		System.out.println("And we are running ...");
		
		optimum = opt.optimize(
				func, GoalType.MINIMIZE, grad, new InitialGuess(rvg.nextVector()), new MaxEval(300));    // runs the optimization    
		
		
		System.out.println("\n\n");
		System.out.println("Number of evaluations: " + opt.getEvaluations());
		System.out.println("Max number of evaluations: " + opt.getMaxEvaluations());
		System.out.println("Number of iterations: " + opt.getIterations());
		System.out.println("Max number of iterations: " + opt.getMaxIterations());
			
	}
	
	public PointValuePair getOptimum() {
		return optimum;
	}
	
}
