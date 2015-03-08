package graph_generators;

import java.util.ArrayList;

import linkpred_batch.RandomWalkGraph;
import linkpred_batch.Ranker;

import org.apache.commons.math3.util.Pair;

import cern.colt.matrix.tdouble.DoubleMatrix1D;

/**
 * 
 * Used for generating artificial graphs, for testing purpose
 *
 */
public abstract class GraphGenerator {
		
	/**
	 * Generate artificial graph used for testing purpose as explained in the paper, 
	 * the graph representation is an ArrayList of FeatureFields 	
	 *  
	 * @param n: number of nodes
	 * @param f: number of features
	 * @param s: the index of the starting node
	 * @param topN: the top ranked N nodes to be put in the D set
	 * @param trueParameters: the parameters used for building the adjacency matrix
	 * @param alpha: the damping factor used for the pagrank when building the D set 
	 * @return RandomWalkGraph
	 */
	public abstract RandomWalkGraph generate (int n, int f, int s);
	
	
	/**
	 * Builds the D set (created links) for synthetic graph and known
	 * parameter values, by taking the first topN highest ranked nodes.
	 * s is the node whose links we are looking at.
	 *
	 * @param graph: the Graph for which the D set is being generated
	 * @param D: the D set
	 * @param L: the L set
	 * @param topN: the top ranked N nodes to be put in the D set
	 * @param parameters: the parameters used for building the adjacency matrix
	 * @param alpha: the damping factor used for the pagerank when building the D set
	 */
	public void buildDandL (RandomWalkGraph graph, int topN, DoubleMatrix1D parameters, double alpha) {
		ArrayList<Pair<Integer, Double>> idRankPairs = Ranker.rankAndSort(
				graph, parameters, alpha);                                                // find pageranks
		
		int i = 0;
		int count = 0;
		while (count < topN) {                                                            // put the highest ranked in D and the rest in L
			if (!graph.hasLink(graph.s, idRankPairs.get(i).getFirst()) &&
					graph.s != idRankPairs.get(i).getFirst()) {				
				graph.D.add(idRankPairs.get(i).getKey());
				count++;
			}
			i++;
		}
		
		while (i < idRankPairs.size()) {
			if (!graph.hasLink(graph.s, idRankPairs.get(i).getFirst())&&
					graph.s != idRankPairs.get(i).getFirst()) 
				graph.L.add(idRankPairs.get(i).getKey());	
			i++;
		}			
	}
	
}
