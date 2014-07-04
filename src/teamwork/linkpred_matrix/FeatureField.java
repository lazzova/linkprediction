package teamwork.linkpred_matrix;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;

public class FeatureField {
	int row;
	int column;
	DoubleMatrix1D features;
	
	/**
	 * Constructor
	 * 
	 * @param row
	 * @param column
	 * @param features
	 */
	public FeatureField(int row, int column, double[] features) {
		super();
		this.row = row;
		this.column = column;
		this.features = new DenseDoubleMatrix1D(features);
	}
	
	/**
	 * Constructor
	 */
	public FeatureField () {}

}
