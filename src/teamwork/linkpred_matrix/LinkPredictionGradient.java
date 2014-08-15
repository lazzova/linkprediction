package teamwork.linkpred_matrix;

import org.apache.commons.math3.analysis.MultivariateVectorFunction;

public class LinkPredictionGradient implements MultivariateVectorFunction {
	private LinkPrediction lp;    // the same object as in LinkPredictionCost
	
	public LinkPredictionGradient (LinkPrediction lp) {
		this.lp = lp;
	}

	@Override
	public double[] value(double[] point) throws IllegalArgumentException {
		return lp.getGradient();
	}

}
