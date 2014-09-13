package linkpred_batch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import org.apache.commons.math3.util.Pair;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;

/**
 * 
 * Class used for calculating and sorting pageranks
 *
 */
public class Ranker {

	/**
	 * Calculates the pagerank and sorts the nodes with respect to it
	 * 
	 * @param graph: the graph
	 * @param parameters: the parameters used by the weighting function for building the adjacency matrix
	 * @param alpha: damping factor
	 * @return ArrayList<Pair<Integer, Double>>
	 */
	public static ArrayList<Pair<Integer, Double>> rankAndSort (
			RandomWalkGraph graph, DoubleMatrix1D parameters, double alpha) {
		graph.buildAdjacencyMatrix(parameters);
		SparseCCDoubleMatrix2D Qt = graph.buildTransitionTranspose(alpha);
		DoubleMatrix1D pagerank = pagerank(Qt);
		
		// sort the ranks in ascending order
		ArrayList<Pair<Integer, Double>> idRankPairs = new ArrayList<Pair<Integer, Double>>();
		for (int i = 0; i < pagerank.size(); i++) 
			if (i != graph.s && graph.A.get(graph.s, i) == 0)                      // the node is not s, and has no links to s previously
				idRankPairs.add(new Pair<Integer, Double> (i, pagerank.get(i)));
						
		Collections.sort(idRankPairs, new Comparator<Pair<Integer, Double>> () {
		@Override
			public int compare(Pair<Integer, Double> o1,
					Pair<Integer, Double> o2) {
				if (o1.getValue() > o2.getValue()) return 1;
				if (o2.getValue() > o1.getValue()) return -1;
				return 0;
			}		
		});
		Collections.reverse(idRankPairs);

		return idRankPairs;
	}
	
	
	/**
	* Calculates the pagerank, given a transition matrix,
	* using the power method
	* 
	* @param Qt: transpose of the transition probability matrix
	* @return DoubleMatrix1D
	*/
	public static DoubleMatrix1D pagerank (SparseCCDoubleMatrix2D Qt) {
		
		int n = Qt.rows();
		DoubleMatrix1D p = new DenseDoubleMatrix1D(n);              // current iteration
		DoubleMatrix1D oldP = new DenseDoubleMatrix1D(n);           // previous iteration
				
		p.assign(1.0 / n);                                          // pagerank initialization 
		
		do {
		
			oldP.assign(p);
			Qt.zMult(oldP, p);
					
			oldP.assign(p, new DoubleDoubleFunction() {
		
				@Override
				public double apply(double arg0, double arg1) {
					return Math.abs(arg0-arg1);
				}
			});
		
		} while (oldP.zSum() > 1E-6);                                // convergence check
		
		return p;
	}
	
	
	/**
	 * Predicts new links by finding the highest ranked nodes given the learned parameters
	 * (Used after the training)
	 * 
	 * @param graph: the graph
	 * @param parameters: the parameters used by the weighting function for building the adjacency matrix
	 * @param alpha: damping factor
	 * @param linksNumber: the number of links to predict
	 * @return ArrayList<Integer> 
	 */
	public static ArrayList<Integer> predictLinks (
			RandomWalkGraph graph, DoubleMatrix1D parameters, double alpha, int linksNumber) {
		ArrayList<Pair<Integer, Double>> highestRanked = rankAndSort(graph, parameters, alpha);
		ArrayList<Integer> links = new ArrayList<Integer> ();
		int count = 0;
		for (int i = 0; count < linksNumber; i++) {
			if (!graph.hasLink(graph.s, highestRanked.get(i).getFirst())) {				
				links.add(highestRanked.get(i).getFirst());
				count++;
			}
		}
		
		return links;
	}

}
