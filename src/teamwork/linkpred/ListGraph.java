package teamwork.linkpred;

import java.util.ArrayList;

public class ListGraph {
	ArrayList<Edge> [] adjList;                                  // adjacency list
	int n;                                                       // number of nodes
	int m;                           	                         // number of features
	
	@SuppressWarnings("unchecked")
	public ListGraph (int n, int m) {
		this.n = n;
		this.m = m;
		this.adjList = (ArrayList<Edge> []) new ArrayList [n];
		for (int i = 0; i < n; i++) {
			this.adjList[i] = new ArrayList<Edge>();
		}
	}	
	
	public void addEdge (Edge e) {
		adjList[e.v1].add(e);
		adjList[e.v2].add(e);
	}
	
	public void addEdge (int v1, int v2) {
		Edge e = new Edge(v1, v2);
		addEdge(e);
	}
	
	public void addEdge (int from, int to, double [] features) {
		Edge e = new Edge(from, to, features);
		addEdge(e);
	}
	
	public void addEdge (int from, int to, double [] features, double weight) {
		Edge e = new Edge(from, to, features, weight);
		addEdge(e);
	}
	
	public double sumWeights (int index) {
		double sum = 0;
		for (int i = 0; i < adjList[index].size(); i++)
			sum += adjList[index].get(i).weight;
		return sum;
	}
	
	public Edge getEdge (int from, int to) {
		/** Retunrs the Edge between from and to */
		Edge e = null;
		for (int i = 0; i < adjList[from].size(); i++) {
			e = adjList[from].get(i);
			if (e.other(from) == to)
				return e;
		}
		return e;
	}
	
	public boolean hasEdge (Edge e) {
		/** Returns true if Edge e exists */
		for (int i = 0; i < adjList[e.v1].size(); i++) 
			if (adjList[e.v1].get(i).other(e.v1) == e.v2)
				return true;
		
		return false;
	}
	
	@Override
	public String toString() {
		/** For debuging purpose */
		StringBuilder s = new StringBuilder();
		for (int i = 0; i < n; i++) {
			s.append(i);
			s.append(": ");
			for (int j = 0; j < adjList[i].size(); j++) {
				s.append(adjList[i].get(j).other(i));
				s.append(" weight:");
				s.append(adjList[i].get(j).weight);
				s.append(" [");
				for (int k = 0; k < adjList[i].get(j).features.length; k++)
					s.append(String.format("%.2f ", adjList[i].get(j).features[k]));
				s.append("]\n");
			}
			s.append("\n");
		}
		return s.toString();
	}
}
