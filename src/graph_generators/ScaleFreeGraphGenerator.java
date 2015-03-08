package graph_generators;

import java.util.ArrayList;
import java.util.Date;

import org.apache.commons.math3.random.GaussianRandomGenerator;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomVectorGenerator;
import org.apache.commons.math3.random.UncorrelatedRandomVectorGenerator;

import linkpred_batch.FeatureField;
import linkpred_batch.Network;
import linkpred_batch.RandomWalkGraph;

public class ScaleFreeGraphGenerator extends GraphGenerator {
	/**Random number generator*/
	private JDKRandomGenerator rand;
		
	public ScaleFreeGraphGenerator() {
		super();
		this.rand = new JDKRandomGenerator();
		this.rand.setSeed(new Date().getTime());		
	}

	public RandomWalkGraph generate (int n, int f, int s) {
		int [] degCumulative = new int [n];	                                        // array for cumulative degree sums
		int [] deg = new int [n];                                                   // array of node degrees
			
		int k;
		int len;
		int randNum;
		             
		RandomVectorGenerator randomVector = 
				new UncorrelatedRandomVectorGenerator(
				f, new GaussianRandomGenerator(rand));
		
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
				
		return g;
	}

}
