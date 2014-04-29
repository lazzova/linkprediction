package teamwork.linkpred_matrix;

import org.apache.commons.math3.analysis.MultivariateFunction;

public class LinkPredictionCost implements MultivariateFunction {
	private LinkPrediction lp;    // the same object as in LinkPredictionGradient
	
	public LinkPredictionCost (LinkPrediction lp) {
		this.lp = lp;
	}

	@Override
	public double value(double[] point) {
		return lp.getCost(point);
	}

}
