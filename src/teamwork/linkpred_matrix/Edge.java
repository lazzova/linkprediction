package teamwork.linkpred_matrix;

public class Edge {
	double weight;
	double [] features;

    public Edge() {
        this.weight = 0;
        this.features = null;
    }

    public Edge( double[] features) {
        this.weight = 0;
        this.features = features;
    }

    public void setFeatures(double[] features) {
        this.features = features;
    }
    
    
        
        

}
