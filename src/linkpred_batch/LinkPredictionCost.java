package linkpred_batch;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.PointValuePair;

public class LinkPredictionCost implements MultivariateFunction {
	private LinkPredictionTrainer lp;                         // the same object as in LinkPredictionGradient
	
	public PointValuePair optimum;
	
	public LinkPredictionCost (LinkPredictionTrainer lp) {
		this.lp = lp;
		optimum = new PointValuePair(null, Double.MAX_VALUE);
	}

	@Override
	public double value(double[] point) {
		//System.out.println("COST"); // TODO
		double d = 0;
		try {
			d = lp.getCost(point);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
				
		System.out.println("Cost: " + d);
		System.out.println("Parameters: ");
		for (int i = 0; i < point.length; i++) 
			System.out.print(point[i] + " ");
		System.out.println();
		System.out.println();
		
		
		if (d < optimum.getValue()) {
			optimum = new PointValuePair(point, d);
		}
		
		return d;
	}

}
