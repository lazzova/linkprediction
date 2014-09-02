package linkpred_batch;

import org.apache.commons.math3.optim.PointValuePair;

public class GradientDescent {
	private static final String NUMBER_OF_ITERATIONS_EXCEDEED = "Maximum number of iteration excedeed.";
	private static final String GRADIENT_CONVERGERNCE = "Gradient convergence: the gradient is smaller than the given treshold";
	private static final String MINIMAL_COST = "The cost is smaller than the given threshold";
	private static final String INVALID = "Something went wrong";
	
	private LinkPredictionTrainer lp;
	private int maxIterations;
	private double gradientTreshold;
	private double costThreshold;
	
	private String stopReason;
	private PointValuePair opt;

	public GradientDescent(LinkPredictionTrainer lp, int maxIterations, double gradientTreshold, double costThreshold) {
		this.lp = lp;
		this.maxIterations = maxIterations;
		this.gradientTreshold = gradientTreshold;
		this.costThreshold = costThreshold;
		this.stopReason = INVALID;
		this.opt = new PointValuePair(null, Double.MAX_VALUE);
	}

	public PointValuePair optimize (double [] initialParameters) throws InterruptedException {
		double [] parameters = initialParameters;
		double [] gradient = lp.getGradient(parameters);
		double cost = lp.getCost(parameters);
		int iteration = 0;
		
		while (true) {
			if (cost < opt.getValue())
				opt = new PointValuePair(parameters, cost);
			
			if (iteration > maxIterations) {
				stopReason = NUMBER_OF_ITERATIONS_EXCEDEED;
				return opt;
			}
			
			if (cost < costThreshold) {
				stopReason = MINIMAL_COST;
				return opt;
			}
			
			if (gradientConverged(gradient)) {
				stopReason = GRADIENT_CONVERGERNCE;
				return opt;
			}
			
			printCurrentIteration(iteration, cost, parameters, gradient);
			updateParameters(parameters, gradient);
			gradient = lp.getGradient(parameters);
			cost = lp.getCost(parameters);
			iteration++;
		}
	}
	
	
	public PointValuePair multiStartOptimize (int restarts, double [] initialParameters) throws InterruptedException {
		PointValuePair multistartOpt = opt;
		String multistartStopReason = INVALID;
		while (--restarts > 0) {
			System.out.println("\n\n**RESTART**\n\n");
			this.stopReason = INVALID;
			this.opt = new PointValuePair(null, Double.MAX_VALUE);
			this.optimize(initialParameters);
			if (opt.getValue() < multistartOpt.getValue()) {
				multistartOpt = opt;
				multistartStopReason = this.stopReason;
			}
		}		
		this.stopReason = multistartStopReason;
		return multistartOpt;
	}
	
		
	private boolean gradientConverged (double [] gradient) {
		for (int i = 0; i < gradient.length; i++)
			if (Math.abs(gradient[i]) >= gradientTreshold)
				return false;
		return true;
	}
	
	
	private void updateParameters (double [] parameters, double [] gradient) {
		for (int i = 0; i < gradient.length; i++)
			parameters[i] -= gradient[i];
	}

	public String getStopReason() {
		return stopReason;
	}
	
	public void printCurrentIteration (int iteration, double cost, double [] parameters, double [] gradient) {
		System.out.println(iteration + 1);
		System.out.println("Cost: " + cost);
		System.out.println("Parameters: ");
		for (int i = 0; i < parameters.length; i++) 
			System.out.print(parameters[i] + " ");		
		System.out.println();
		for (int i = 0; i < gradient.length; i++) 
			System.out.print(gradient[i] + " ");		
		System.out.println();
		System.out.println();
	}
}
