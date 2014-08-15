package linkpred_stochastic;

import org.apache.commons.math3.analysis.MultivariateFunction;

public class LinkPredictionCostStochastic implements MultivariateFunction {
	private LinkPredictionStochastic lp;                         // the same object as in LinkPredictionGradient
	
	public LinkPredictionCostStochastic (LinkPredictionStochastic lp) {
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
