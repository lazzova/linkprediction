package teamwork.linkpred;

public class Edge {
	double weight;
	double [] features;
	int v1; 
	int v2; 
	
	public Edge (int v1, int v2) {
		this.weight = 0;
		this.features = null;
		this.v1 = v1;
		this.v2 = v2;
	}
	
	public Edge (int v1, int v2, double [] features) {
		this.weight = 0;
		this.features = features;
		this.v1 = v1;
		this.v2 = v2;
	}
	
	public Edge (int v1, int v2, double [] features, double weight) {
		this.weight = weight;
		this.features = features;
		this.v1 = v1;
		this.v2 = v2;
	}
	
	public void setFeatures(double[] features) {
		this.features = features;
	}
	
	public int other (int node) {
		/** return the other node */
		if (node == v1)
			return v2;
		return v1;
	}
}
