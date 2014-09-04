package linkpred_batch;

import org.apache.commons.math3.optim.PointValuePair;

public class GradientDescent {
	/**Stopping reason 1*/
	private static final String NUMBER_OF_ITERATIONS_EXCEDEED = "Maximum number of iteration excedeed.";
	/**Stopping reason 2*/
	private static final String GRADIENT_CONVERGERNCE = "Gradient convergence: the gradient is smaller than the given treshold";
	/**Stopping reason 3*/
	private static final String MINIMAL_COST = "The cost is smaller than the given threshold";
	/**Stopping reason 4*/
	private static final String INVALID = "Something went wrong";
	
	/**The trainer object*/
	private LinkPredictionTrainer lp;
	/**Maximum number of iterations allowed*/
	private int maxIterations;
	/**Minimum gradient threshold*/
	private double gradientThreshold;
	/**Minimum cost*/
	private double costThreshold;
	
	
	/**Stopping reason*/
	private String stopReason;
	/**Optimal point and value*/
	private PointValuePair opt;

	
	/**
	 * Constructor
	 * 
	 * @param lp
	 * @param maxIterations
	 * @param gradientThreshold
	 * @param costThreshold
	 */
	public GradientDescent(LinkPredictionTrainer lp, int maxIterations, double gradientThreshold, double costThreshold) {
		this.lp = lp;
		this.maxIterations = maxIterations;
		this.gradientThreshold = gradientThreshold;
		this.costThreshold = costThreshold;
		this.stopReason = INVALID;
		this.opt = new PointValuePair(null, Double.MAX_VALUE);
	}

	
	/**
	 * Gradient descent based optimization function
	 * 
	 * @param initialParameters
	 * @return PointValuePair
	 * @throws InterruptedException
	 */
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
	
	
	/**
	 * Gradient descent multistart optimization function
	 * 
	 * @param restarts
	 * @param initialParameters
	 * @return PointValuePair
	 * @throws InterruptedException
	 */
	public PointValuePair multiStartOptimize (int restarts, double [] initialParameters) throws InterruptedException {
		PointValuePair multistartOpt = opt;
		String multistartStopReason = INVALID;
		while (--restarts > 0) {
			System.out.println("\n\n**RESTART**\n\n");
			this.stopReason = INVALID;
			this.opt = new PointValuePair(null, Double.MAX_VALUE);
			for (int i = 0; i < initialParameters.length; i++)
				initialParameters[i] = Math.random() * 2 - 1;  // TODO random within interval [-1,1];
			
			this.optimize(initialParameters);
			if (opt.getValue() < multistartOpt.getValue()) {
				multistartOpt = opt;
				multistartStopReason = this.stopReason;
			}
			
			if (multistartStopReason.equals(MINIMAL_COST))
				break;
		}		
		this.stopReason = multistartStopReason;
		return multistartOpt;
	}
	
	
	/**
	 * Checks if the gradient has reached the threshold
	 * 
	 * @param gradient
	 * @return boolean
	 */
	private boolean gradientConverged (double [] gradient) {
		for (int i = 0; i < gradient.length; i++)
			if (Math.abs(gradient[i]) >= gradientThreshold)
				return false;
		return true;
	}
	
	
	/**
	 * Updates the parameters 
	 * 
	 * @param parameters
	 * @param gradient
	 */
	private void updateParameters (double [] parameters, double [] gradient) {
		for (int i = 0; i < gradient.length; i++)
			parameters[i] -= gradient[i];
	}

	
	/**
	 * Returns the stopping reason
	 * 
	 * @return String
	 */
	public String getStopReason() {
		return stopReason;
	}
	
	
	/**
	 * Prints info regarding the current iteration (for debugging purpose)
	 * 
	 * @param iteration
	 * @param cost
	 * @param parameters
	 * @param gradient
	 */
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
