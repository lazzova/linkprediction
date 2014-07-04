package teamwork.linkpred_matrix;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ml.distance.ManhattanDistance;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.util.Pair;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.function.tdouble.DoubleFunction;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.solver.DefaultDoubleIterationMonitor;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;

public class Main {
	
	/**
	 * Main
	 *  
	 * @param args
	 */
	public static void main(String[] args) {
		
		long start = System.nanoTime();
		System.out.println("Graph generation start");            //TODO
		
		int g = 50;                                              // number of graphs
		int n = 10000;                                           // number of nodes
		int f = 2;                                               // number of features
		
		GraphGeneration.initRandom(f);                                     // build the graph
		Graph [] featureMat = new Graph [g];
		for (int i = 0; i < g; i++)
			featureMat[i] = GraphGeneration.generate(n);
		
		System.out.println("Graph generation end");				 //TODO
		
		int s = 0;                                               // the node whose links we learn
		double alpha = 0.2;                                      // damping factor
		double b = 1; //1e-6;                                    // WMW function parameter
		double lambda = 1;                                       // regularization parameter 
		double [] param = {1, -1};                               // parameters vector
		DoubleMatrix1D parameters = new DenseDoubleMatrix1D(param);
		
		
		System.out.println("Building D start");                  //TODO
		
		int topRanked = 15;                                      // building D set (linked set)
		for (int i = 0; i < g; i++)
			featureMat[i].buildD(topRanked, parameters, s, alpha);
				
		//TODO
		long end = System.nanoTime();
		System.out.println("Building D end");
		System.out.println("Graph generation finished in " + (end-start)/1E9 + " seconds");
		System.out.println("Memory used: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1E6 + "MB");
				
		
		/*
		LinkpredProblem problem = new LinkpredProblem(graphs, connected, s, D, alpha, lambda, b);
		problem.optimize();
		PointValuePair optimum = problem.getOptimum();
		
		System.out.println("Function minimum: " + optimum.getValue() + "\nParameters: " + 
		        optimum.getPoint()[0] + " " + optimum.getPoint()[1]);
		*/
	}
	
	
	/*
	/**
	 * Builds the D set (created links) for synthetic graph and known
	 *  parameter values, by taking the first topN highest ranked nodes.
	 *  s is the node whose links we are looking at
	 * 
	 * @param topN 
	 * @param graph
	 * @param trueParameters
	 * @param s
	 * @param alpha
	 * @param D will contain the linked set, assumed empty at start
	 * @param L will contain the no-link set, assumed empty at start
	 * @return
	 *
	public static void buildD (int topN, FeatureMatrix graph, DoubleMatrix1D trueParameters, 
								  int s, double alpha, ArrayList<Pair<Integer, Double>> D, 
			                      ArrayList<Pair<Integer, Double>> L) {
		
		// array of indices
		int n = graph.dim;
		int [] nodes = new int [n];
		for (int i = 1; i < n; i++)
			nodes[i] = i;
		
		// find pageranks
		SparseCCDoubleMatrix2D A = graph.toAdjacencyMatrix(trueParameters);
		DoubleMatrix2D Q = buildTransitionMatrix(A, graph, s, alpha);
		DoubleMatrix1D rank = pagerank(Q);
		
		// sort the ranks in ascending order
		for (int i = 0; i < rank.size(); i++)
			L.add(new Pair<Integer, Double> (i, rank.get(i)));
		Collections.sort(L, new Comparator<Pair<Integer, Double>> () {

			@Override
			public int compare(Pair<Integer, Double> o1,
					Pair<Integer, Double> o2) {
				if (o1.getValue() > o2.getValue()) return 1;
				if (o2.getValue() > o1.getValue()) return -1;
				return 0;
			}
			
		});
		
		// put the highest ranked in D and remove those from L
		Pair<Integer, Double> pair;
		for (int i = 0; D.size() < topN; i++) {
			pair = L.get(L.size()-i-1);
			if (pair.getKey() != s && A.get(s, pair.getKey()) == 0) {      // the node is not s, and has no links to s previously
				D.add(pair);
				L.remove(L.size()-i-1);
			}
		}
		
			
	}
	
	
	/**
	 * Builds the transition matrix for given adjacency matrix
	 * and s as starting node
	 *
	 * @param A
	 * @param s
	 * @param alpha
	 * @return DoubleMatrix2D
	 *
	private static DoubleMatrix2D buildTransitionMatrix (
	          SparseCCDoubleMatrix2D A, FeatureMatrix fm, int s, double alpha) {

		SparseCCDoubleMatrix2D Q = new SparseCCDoubleMatrix2D(A.rows(), A.columns());
	
		// (1-alpha) * A[i][j] / sumElements(A[i])) + 1(j == s)
		double value;
		int r, c;
		for (int i = 0; i < fm.list.size(); i++) {
			r = fm.list.get(i).row;
			c = fm.list.get(i).column;
			value = A.get(r, c);
			value *= (1 - alpha);
			value /= A.viewRow(r).zSum();
			Q.set(r, c, value);
		
			if (r == c) continue;
		
			value = A.get(c, r);
			value *= (1 - alpha);
			value /= A.viewRow(c).zSum();
			Q.set(c, r, value);
		}
				
		for (int i = 0; i < A.rows(); i++) {
			value = A.get(i, s);
			value += alpha;
			Q.set(i, s, value);
		}
	
		return Q;
	}

	/**
	 * Calculates the pagerank, given a transition matrix,
	 * using the power method
	 * 
	 * @param Q
	 * @return
	 *
	private static DoubleMatrix1D pagerank (DoubleMatrix2D Q) {
		
		int n = Q.rows();
		DoubleMatrix1D p = new DenseDoubleMatrix1D(n);           // current iteration
		DoubleMatrix1D oldP = new DenseDoubleMatrix1D(n);        // previous iteration
		DenseDoubleAlgebra algebra = new DenseDoubleAlgebra();
		DoubleMatrix2D Qtranspose = algebra.transpose(Q);
		
		p.assign(1.0 / n);                                       // pagerank initialization 
		
		do {
			
			for (int i = 0; i < n; i++)
				oldP.set(i, p.get(i));
			p = algebra.mult(Qtranspose, p);
			
			oldP.assign(p, new DoubleDoubleFunction() {
				
				@Override
				public double apply(double arg0, double arg1) {
					return Math.abs(arg0-arg1);
				}
			});
			
		} while (algebra.norm1(oldP) > 1E-6);                    // convergence check
		
		return p;
	}*/
}




