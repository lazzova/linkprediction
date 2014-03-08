package teamwork.linkpred;

public class Edge {
	double weight;
	double [] features;
	int from; // start node
	int to;   // end node
	
	public Edge (int from, int to) {
		this.weight = 0;
		this.features = null;
		this.from = from;
		this.to = to;
	}
	
	public Edge (int from, int to, double [] features) {
		this.weight = 0;
		this.features = features;
		this.from = from;
		this.to = to;
	}
	
	public void setFeatures(double[] features) {
		this.features = features;
	}	
}
