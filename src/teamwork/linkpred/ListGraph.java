package teamwork.linkpred;

import java.util.ArrayList;

public class ListGraph {
	private ArrayList<Edge> [] adjList;   // adjacency list
	private int n;                        // number of nodes
		
	
	@SuppressWarnings("unchecked")
	public ListGraph (int n) {
		this.n = n;
		this.adjList = (ArrayList<Edge> []) new ArrayList [n];
		for (int i = 0; i < n; i++) {
			this.adjList[i] = new ArrayList<Edge>();
		}
	}	
	
	public ArrayList<Edge>[] getAdjList() {
		return adjList;
	}

	public int getN() {
		return n;
	}
	
	public void addEdge (Edge e) {
		adjList[e.getFrom()].add(e);
		adjList[e.getTo()].add(e);
	}
	
	public void addEdge (int from, int to) {
		Edge e = new Edge(from, to);
		addEdge(e);
	}
	
	public void addEdge (int from, int to, byte weightFunction) {
		Edge e = new Edge(from, to, weightFunction);
		addEdge(e);
	}
	
	public void addEdge (int from, int to, byte weightFunction, double [] features) {
		Edge e = new Edge(from, to, weightFunction, features);
		addEdge(e);
	}
	
	public double sumWeights (int index) {
		double sum = 0;
		for (int i = 0; i < adjList[index].size(); i++)
			sum += adjList[index].get(i).getWeight();
		return sum;
	}
}
