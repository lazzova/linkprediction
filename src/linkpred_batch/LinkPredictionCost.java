package linkpred_batch;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.PointValuePair;

public class LinkPredictionCost implements MultivariateFunction {
	/**The training object, the same as in LinkPredictionGradient*/
	private LinkPredictionTrainer lp;                         
	/**Optimal point and value*/
	public PointValuePair optimum;
	
	
	/**
	 * Constructor
	 * 
	 * @param lp: the trainer object
	 */
	public LinkPredictionCost (LinkPredictionTrainer lp) {
		this.lp = lp;
	}

	
	/**
	 *  Returns the value of the cost function,
	 *  given parameters
	 *  
	 *  @param point: the parameters
	 */
	@Override
	public double value(double[] point) {
		double d = lp.getCost(point);
						
		System.out.println("Cost: " + d);
		System.out.println("Parameters: ");
		for (int i = 0; i < point.length; i++) 
			System.out.print(point[i] + " ");
		System.out.println();
		System.out.println();
				
		return d;
	}

}
