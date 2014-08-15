package linkpred_stochastic;

import org.apache.commons.math3.analysis.MultivariateVectorFunction;

public class LinkPredictionGradientStochastic implements MultivariateVectorFunction {
	private LinkPredictionStochastic lp;               // the same object as in LinkPredictionCost
	int gCounter;                                      // current graph
	int g;                                             // number of graphs
	
	public LinkPredictionGradientStochastic (LinkPredictionStochastic lp) {
		this.lp = lp;
		this.gCounter = -1;
		this.g = lp.g;
	}

	@Override
	public double[] value(double[] point) throws IllegalArgumentException {
		gCounter++;
		if (gCounter >= g)
			gCounter = 0;
		
		double[] grad = null;
		try {
			grad = lp.getGradient(point, gCounter);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		return grad;
	}

}
