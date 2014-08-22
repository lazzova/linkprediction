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
		try {
			grad = lp.getGradient(point);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		
		System.out.println("Point: ");
		System.out.print(point[0] + "\t" + point[1]);
		System.out.println();
		System.out.println("Gradient: ");
		System.out.print(grad[0] + "\t" + grad[1]);
		System.out.println("\n");
		
		boolean tooSmall = true;
		for (int i = 0; i < grad.length; i++) {
			if (grad[i] > 0.001) {
				tooSmall = false;
				break;
			}
		}
			
		if (tooSmall) grad = new double[point.length];
		
		return grad;
	}

}
