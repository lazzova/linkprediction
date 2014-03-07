package teamwork.linkpred;

import java.util.HashMap;

public class LinkPrediction {
	private ListGraph graph;
	private double alpha;               // damping factor
	private double lambda;              // regularization parameter
	private double [] p;                // page rank
	private double [] dp;               // page rank gradient
	private double J;                   // cost
	private double gradient;            // gradient of J
	private double wPredicted;          // predicted parameters
	private HashMap<String, Double> q;  // transition matrix
	// We call it the transition matrix 
    // by convention but we represent it 
	// with a HashMap
	private double [] dTransition;       // transition gradient TODO
	private int s;                       // index of the s node
	
	public LinkPrediction () {
	    q = new HashMap<String,Double> ();         		
	}

	public ListGraph getGraph() {
		return graph;
	}

	public void setGraph(ListGraph graph) {
		this.graph = graph;
	}

	public double getAlpha() {
		return alpha;
	}

	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}

	public double getLambda() {
		return lambda;
	}

	public void setLambda(double lambda) {
		this.lambda = lambda;
	}

	public double[] getP() {
		return p;
	}

	public void setP(double[] p) {
		this.p = p;
	}

	public double[] getDp() {
		return dp;
	}

	public void setDp(double[] dp) {
		this.dp = dp;
	}

	public double getJ() {
		return J;
	}

	public void setJ(double j) {
		J = j;
	}

	public double getGradient() {
		return gradient;
	}

	public void setGradient(double gradient) {
		this.gradient = gradient;
	}

	public double getWPredicted() {
		return wPredicted;
	}

	public void setWPredicted(double w_predicted) {
		this.wPredicted = w_predicted;
	}

	public double getwPredicted() {
		return wPredicted;
	}

	public void setwPredicted(double wPredicted) {
		this.wPredicted = wPredicted;
	}

	public HashMap<String, Double> getQ() {
		return q;
	}
	
	public void buildTransitionMatrix () {
		/** Builds the transition matrix q */
		for (int i = 0; i < graph.getN(); i++) {
			for (int j = 0; j < graph.getAdjList()[i].size(); i++) {
				q.put(String.format("%d,%d", i, j),
					  graph.getAdjList()[i].get(j).getWeight() / 
					  graph.sumWeights(i));
			}
		}
	}

		
}
