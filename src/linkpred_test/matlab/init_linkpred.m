function [obj] = init_linkpred (g, n, fs, s, alpha, b, lambda, param, learningrate, interlayer, topn)
    % import all java classes needed    
    import teamwork.linkpred_matrix.*;
    import cern.colt.matrix.tdouble.DoubleMatrix1D;
    import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
    import cern.colt.function.tdouble.DoubleDoubleFunction;
    import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;
    import cern.jet.math.tdouble.DoubleFunctions;
    import org.apache.commons.math3.util.Pair;
    import org.apache.commons.math3.random.GaussianRandomGenerator;
    import org.apache.commons.math3.random.JDKRandomGenerator;
    import org.apache.commons.math3.random.RandomVectorGenerator;
    import org.apache.commons.math3.random.UncorrelatedRandomVectorGenerator;
    import edu.emory.mathcs.csparsej.tdouble.Dcs_util;
    import edu.emory.mathcs.utils.ConcurrencyUtils;
    
    % create an object and initialize the linkprediction problem
    obj = MatlabOptFunction;
    obj.initProblem(g, n, f, s, alpha, b, lambda, param, learningrate, interlayer, topn);        
end
