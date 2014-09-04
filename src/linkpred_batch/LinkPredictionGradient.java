package linkpred_batch;

import org.apache.commons.math3.analysis.MultivariateVectorFunction;


public class LinkPredictionGradient implements MultivariateVectorFunction {
	/**The training object, the same as in LinkPredictionCost*/
	private LinkPredictionTrainer lp;    
			
	/**
	 * Constructor
	 * 
	 * @param lp: the trainer object
	 */
	public LinkPredictionGradient (LinkPredictionTrainer lp) {
		this.lp = lp;		
	}

	
	/**
	 *  Returns the gradient of the cost function,
	 *  given parameters
	 *  
	 *  @param point: the parameters
	 */
	@Override
	public double[] value(double[] point) throws IllegalArgumentException {
		
		double[] grad = null;
		try {
			grad = lp.getGradient(point);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		
		System.out.println("Point: ");
		for (int i = 0; i < point.length; i++)
			System.out.print(point[i] + "\t");
		System.out.println();
		System.out.println("Gradient: ");
		for (int i = 0; i < grad.length; i++)
			System.out.print(grad[i] + "\t");
		System.out.println("\n");
		
		/*
		boolean tooSmall = true;
		for (int i = 0; i < grad.length; i++) {
			if (grad[i] > 0.00001) {
				tooSmall = false;
				break;
			}
		}
		
			
		if (tooSmall || lp.getCost(point) < 2) grad = new double[point.length];
		*/
		
		return grad;
	}

}
