package linkpred_batch;

import java.util.ArrayList;
import java.util.Date;

import org.apache.commons.math3.random.GaussianRandomGenerator;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomVectorGenerator;
import org.apache.commons.math3.random.UncorrelatedRandomVectorGenerator;
import org.apache.commons.math3.util.Pair;

import cern.colt.matrix.tdouble.DoubleMatrix1D;


public class ArtificialGraphGenerator {
	/**Random number generator*/
	public static JDKRandomGenerator rand = new JDKRandomGenerator();
	/**Random vector generator*/
	public static RandomVectorGenerator randomVector = null;
	
	
	/**
	 * Set seed and initialize the generator of
	 * random vectors of length f 
	 *
	 * @param f: the number of features
	 */
	public static void initialize (int f) {
		rand.setSeed(new Date().getTime());
		randomVector = 
				new UncorrelatedRandomVectorGenerator(
				f, new GaussianRandomGenerator(rand));
	}
	
	
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
	public static RandomWalkGraph generate (int n, int f, int s, int topN, 
			DoubleMatrix1D trueParameters, double alpha) {
		int [] degCumulative = new int [n];	                                        // array for cumulative degree sums
		int [] deg = new int [n];                                                   // array of node degrees
			
		int k;
		int len;
		int randNum;
		             
		
		ArrayList<FeatureField> featureList = new ArrayList<FeatureField>();
		featureList.add(new FeatureField (0, 1, randomVector.nextVector()));        // connect first three nodes in a triad
		featureList.add(new FeatureField (1, 2, randomVector.nextVector()));
		featureList.add(new FeatureField (2, 0, randomVector.nextVector()));
		deg[0] = 2;
		deg[1] = 2;
		deg[2] = 2;
					
		FeatureField ff;
		for (int i = 3; i < n; i++) {
			for (int j = 0; j < 3; j++) {                                           // generate three links
				if (rand.nextInt(11) < 8) {                                         // select destination node randomly
					randNum = rand.nextInt(i);
					ff = new FeatureField (i, randNum, randomVector.nextVector());
					if (!featureList.contains(ff))
						featureList.add(ff);
					deg[i]++;
					deg[randNum]++;
				}
				
				else {                                                              // select destination node proportionally to its degree
					degCumulative[0] = deg[0];
					for (k = 1; k < i; k++)
						degCumulative[k] = degCumulative[k-1] + deg[k];
					len = k-1;
				    randNum = rand.nextInt(degCumulative[len-1] + 1);	
					k = 0;
				    while (randNum > degCumulative[k])
				    	k++;
				    
				    ff = new FeatureField (i, k, randomVector.nextVector());
				    if (!featureList.contains(ff))
						featureList.add(ff);
				    deg[i]++;
				    deg[k]++;					    
				}								
			}			
		}	
		
		RandomWalkGraph g = new Network(n, f, s, featureList, 
				new ArrayList<Integer>(), new ArrayList<Integer>());
		buildDandL(g, topN, trueParameters, alpha);                                 // build D and L set
		
		return g;
	}
	
	
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
	public static void buildDandL (RandomWalkGraph graph, int topN, DoubleMatrix1D parameters, double alpha) {
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
