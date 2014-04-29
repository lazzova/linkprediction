package teamwork.linkpred_matrix;

import java.util.Date;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
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

public class Main {

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
		*/
	}
}


class Function implements MultivariateFunction {

	@Override
	public double value(double[] param) {
		double x1 = param[0];
		double x2 = param[1];
		return x1*x1 + 2*x2*x2 - 3*x1 + 7*x2 -1;
	}
	
}

class Gradient implements MultivariateVectorFunction {

	@Override
	public double[] value(double[] param) throws IllegalArgumentException {
		return new double [] {2*param[0]-3, 4*param[1]+7};
	}
	
}



