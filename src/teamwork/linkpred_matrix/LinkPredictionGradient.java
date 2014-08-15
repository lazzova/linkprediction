package teamwork.linkpred_matrix;

import org.apache.commons.math3.analysis.MultivariateVectorFunction;

public class LinkPredictionGradient implements MultivariateVectorFunction {
	private LinkPrediction lp;    // the same object as in LinkPredictionCost
	
	public LinkPredictionGradient (LinkPrediction lp) {
		this.lp = lp;
	}

	@Override
	public double[] value(double[] point) throws IllegalArgumentException {
		//System.out.println("gradient");
		double [] grad = lp.getGradient();
		//System.out.println();
		//for (int i = 0; i < grad.length; i++) 
		//	System.out.print(grad[i] + " ");
		//System.out.println();
		return grad;
	}

}
