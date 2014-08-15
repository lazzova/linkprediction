package linkpred_batch;

import org.apache.commons.math3.analysis.MultivariateFunction;

public class LinkPredictionCost implements MultivariateFunction {
	private LinkPrediction lp;                         // the same object as in LinkPredictionGradient
	
	public LinkPredictionCost (LinkPrediction lp) {
		this.lp = lp;
	}

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
