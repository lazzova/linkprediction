package linkpred_batch;

import org.apache.commons.math3.analysis.MultivariateVectorFunction;


public class LinkPredictionGradient implements MultivariateVectorFunction {
	private LinkPredictionTrainer lp;    // the same object as in LinkPredictionCost
		
	public LinkPredictionGradient (LinkPredictionTrainer lp) {
		this.lp = lp;
	}

	@Override
	public double[] value(double[] point) throws IllegalArgumentException {
		
		double[] grad = null;
		//try {
			grad = lp.getGradient(point);
		//} catch (InterruptedException e) {
		//	e.printStackTrace();
		//}
		return grad;
	}

}
