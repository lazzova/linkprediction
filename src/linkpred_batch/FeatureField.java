package linkpred_batch;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;

public class FeatureField {
	int row;                                 // row index of the matrix field
	int column;                              // column index of the matrix field
	DoubleMatrix1D features;                 // the features of the field
	
	
	/**
	 * Constructor
	 * 
	 * @param row: row index of the matrix field
	 * @param column: column index of the matrix field
	 * @param features: the features of the field
	 */
	public FeatureField(int row, int column, double[] features) {
		super();
		this.row = row;
		this.column = column;
		this.features = new DenseDoubleMatrix1D(features);
	}
	
	
	/**
	 * Constructor
	 * 
	 * @param row: row index of the matrix field
	 * @param column: column index of the matrix field
	 * @param features: the features of the field
	 */
	public FeatureField(int row, int column, DoubleMatrix1D features) {
		super();
		this.row = row;
		this.column = column;
		this.features = features;
	}
	
	
	/**
	 * Constructor
	 */
	public FeatureField () {}
	
	
	@Override
	public boolean equals(Object obj) {
		if (obj == null) return false;
		if (this.getClass() != obj.getClass()) return false;
		return column == ((FeatureField) obj).column && 
				row == ((FeatureField) obj).row;
	}
}
