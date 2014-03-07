package teamwork.linkpred;

public class Edge {
	private static byte weightFunction;
	
	// edge-strength functions
	private static final byte EXPONENTIAL = 1; //default
	private static final byte LOGISTIC = 2;
	
	private int id;
	private double weight;
	private double [] features;
	private int from; // start node
	private int to;   // end node
	
	public Edge (int from, int to) {
		this.weight = 0;
		if (Edge.weightFunction == 0)
			Edge.weightFunction = EXPONENTIAL;
		this.features = null;
		this.from = from;
		this.to = to;
	}
	
	public Edge (int from, int to, byte weightFunction) {
		this.weight = 0;
		if (Edge.weightFunction == 0)
			Edge.weightFunction = weightFunction;
		this.features = null;
		this.from = from;
		this.to = to;
	}
	
	public Edge (int from, int to, byte weightFunction, double [] features) {
		this.weight = 0;
		if (Edge.weightFunction == 0)
			Edge.weightFunction = weightFunction;
		this.features = features;
		this.from = from;
		this.to = to;
	}
	
	public int getID () {
		return this.id;
	}
	
	public double getWeight () {
		return this.weight;	
	}
	
	public int getFrom () {
		return this.from;
	}
	
	public int getTo () {
		return this.to;
	}
	
	public double[] getFeatures() {
		return features;
	}

	public void setFeatures(double[] features) {
		this.features = features;
	}

	public double setWeigth (double [] parameters) {
		if (this.weightFunction == EXPONENTIAL)
			return exponental (dotProduct(parameters));
		else
			return sigmoid (dotProduct(parameters));
	}
	
	private double dotProduct (double [] parameters) {
		double dp = 0;
		for (int i = 0; i < this.features.length; i++)
			dp += (this.features[i] * parameters[i]);
		return dp;
	}
	
	private double exponental (double z) {
		return Math.exp(z);
	}
	
	private double sigmoid (double z) {
		return 1.0 / (1+ Math.exp(-z));
	}
}
