package teamwork.linkpred;

import java.util.ArrayList;

public class ListGraph {
	ArrayList<Edge> [] adjList;   // adjacency list
	int n;                        // number of nodes
		
	
	@SuppressWarnings("unchecked")
	public ListGraph (int n) {
		this.n = n;
		this.adjList = (ArrayList<Edge> []) new ArrayList [n];
		for (int i = 0; i < n; i++) {
			this.adjList[i] = new ArrayList<Edge>();
		}
	}	
	
	public void addEdge (Edge e) {
		adjList[e.from].add(e);
		adjList[e.to].add(e);
	}
	
	public void addEdge (int from, int to) {
		Edge e = new Edge(from, to);
		addEdge(e);
	}
	
	public void addEdge (int from, int to, double [] features) {
		Edge e = new Edge(from, to, features);
		addEdge(e);
	}
	
	public double sumWeights (int index) {
		double sum = 0;
		for (int i = 0; i < adjList[index].size(); i++)
			sum += adjList[index].get(i).weight;
		return sum;
	}
}
