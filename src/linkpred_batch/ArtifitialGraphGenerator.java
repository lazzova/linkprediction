package linkpred_batch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;

import org.apache.commons.math3.random.GaussianRandomGenerator;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomVectorGenerator;
import org.apache.commons.math3.random.UncorrelatedRandomVectorGenerator;
import org.apache.commons.math3.util.Pair;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;

public class ArtifitialGraphGenerator {
	public static JDKRandomGenerator rand = new JDKRandomGenerator();
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
	 * Generate graph in FeatureMatrix model
	 * Used for testing purpose	
	 *  
	 * @param n: number of nodes
	 * @param f: number of features
	 * @param s: the index of the starting node
	 * @param topN: the top ranked N nodes to be put in the D set
	 * @param trueParameters: the parameters used for building the adjacency matrix
	 * @param alpha: the damping factor used for the pagrank when building the D set 
	 * @return Graph
	 */
	public static RandomWalkGraph generate (int n, int f, int s, int topN, DoubleMatrix1D trueParameters, double alpha) {
		int [] degCumulative = new int [n];	                     // array for cumulative degree sums
		int [] deg = new int [n];                                // array of node degrees
			
		int k;
		int len;
		int randNum;
		             
		
		ArrayList<FeatureField> featureList = new ArrayList<FeatureField>();
		featureList.add(new FeatureField (0, 1, randomVector.nextVector()));          // connect first three nodes in a triad
		featureList.add(new FeatureField (1, 2, randomVector.nextVector()));
		featureList.add(new FeatureField (2, 0, randomVector.nextVector()));
		deg[0] = 2;
		deg[1] = 2;
		deg[2] = 2;
					
		FeatureField ff;
		for (int i = 3; i < n; i++) {
			for (int j = 0; j < 3; j++) {                        // generate three links
				if (rand.nextInt(11) < 8) {                      // select destination node randomly
					randNum = rand.nextInt(i);
					ff = new FeatureField (i, randNum, randomVector.nextVector());
					if (!featureList.contains(ff))
						featureList.add(ff);
					deg[i]++;
					deg[randNum]++;
				}
				
				else {                                           // select destination node proportionally to its degree
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
		
		RandomWalkGraph g = new Network(n, f, s, featureList, new ArrayList<Integer>(), new ArrayList<Integer>());
		/*g.dim = n;
		g.f = f;
		g.s = s;
		g.list = featureList;
		g.D = new ArrayList<Integer>();
		g.L = new ArrayList<Integer> ();
		g.A = new SparseCCDoubleMatrix2D(n, n);*/
		
		buildDandL(g, topN, trueParameters, alpha);
		
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
	 * @param alpha: the damping factor used for the pagrank when building the D set
	 */
	public static void buildDandL (RandomWalkGraph graph, int topN, DoubleMatrix1D parameters, double alpha) {
		// find pageranks
		ArrayList<Pair<Integer, Double>> idRankPairs = Ranker.predict(graph, parameters, alpha);
		
		// put the highest ranked in D and remove those from L
		int i = 0;
		while (i < topN)
			graph.D.add(idRankPairs.get(i++).getKey());
		while (i < idRankPairs.size())
			graph.L.add(idRankPairs.get(i++).getKey());			
	}		
}
