package linkpred_batch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import org.apache.commons.math3.util.Pair;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;


public class Ranker {

	public static ArrayList<Pair<Integer, Double>> predict (RandomWalkGraph graph, DoubleMatrix1D parameters, double alpha) {
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
	* @return
	*/
	public static DoubleMatrix1D pagerank (SparseCCDoubleMatrix2D Qt) {
		
		int n = Qt.rows();
		DoubleMatrix1D p = new DenseDoubleMatrix1D(n);           // current iteration
		DoubleMatrix1D oldP = new DenseDoubleMatrix1D(n);        // previous iteration
				
		p.assign(1.0 / n);                                       // pagerank initialization 
		
		do {
		
			oldP.assign(p);
			Qt.zMult(oldP, p);
					
			oldP.assign(p, new DoubleDoubleFunction() {
		
				@Override
				public double apply(double arg0, double arg1) {
					return Math.abs(arg0-arg1);
				}
			});
		
		} while (oldP.zSum() > 1E-6);                    // convergence check
		
		return p;
	}

}
