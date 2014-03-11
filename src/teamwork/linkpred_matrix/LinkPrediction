import java.util.HashMap;

/**
 *
 * @author blagoja
 */
public class LinkPrediction {
    static final byte EXPONENTIAL = 1; //default
	static final byte LOGISTIC = 2;
    	
	private MatrixGraph graph;
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
	private byte weightFunction; 
        
        public LinkPrediction () {
	    q = new HashMap<String,Double> ();         		
	}
        
        public void setWeightFunction (byte weightFunction) {
		/** Sets the edge-strength function */
		this.weightFunction = weightFunction;
	}
        
        public void setEdgeWeigths (double [] parameters) {
		/** Sets the weights for all edges in the graph according
		 *  to the chosen edge strength function
		 */
		for (int i = 0; i < graph.n; i++) {
			for (int j = 0; j < graph.n; j++) {
				if (this.weightFunction == EXPONENTIAL)
					graph.adjMatrix[i][j].weight = exponental (dotProduct(
							graph.adjMatrix[i][j].features, parameters));
				else if (this.weightFunction == LOGISTIC)
					graph.adjMatrix[i][j].weight = logistic (dotProduct(
							graph.adjMatrix[i][j].features, parameters));
			}
		}
	}

	public static double dotProduct (double [] features, double [] parameters) {
		/** Calculate dot product between given features vector
		 *  and given parameters vector
		 */
		double dProd = 0;
		for (int i = 0; i < features.length; i++)
			dProd += (features[i] * parameters[i]);
		return dProd;
	}
        
        private double exponental (double z) {
		/** Calculates the exponential function */
		return Math.exp(z);
	}

	private double logistic (double z) {
		/** Calculates the logistic function */
		return 1.0 / (1+ Math.exp(-z));
	}

	public void setAlpha(double alpha) {
		/** Sets the damping factor */
		this.alpha = alpha;
	}
        public void buildTransitionMatrix () {
		/** Builds the transition matrix q */
		double tmp = 0;
		for (int i = 0; i < graph.n; i++) {
			for (int j = 0; j < graph.n; i++) {
				tmp = (1-alpha) * graph.adjMatrix[i][j].weight / 
						  graph.sumWeights(i);
				if (j == s)
					tmp += alpha;
				q.put(String.format("%d,%d", i, j), tmp);
			}
		}
	}	

	public void calculatePageRank () {
		/** Calculates the page rank, using pagerank
		 *  with restarts algorithm
		 */
	}

	public MatrixGraph getGraph() {
		return graph;
	}

	public void setGraph(MatrixGraph graph) {
		this.graph = graph;
	}

	public double getAlpha() {
		return alpha;
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


	public void setwPredicted(double wPredicted) {
		this.wPredicted = wPredicted;
	}

	public HashMap<String, Double> getQ() {
		return q;
	}
}
