public class VertexStack implements java.io.Serializable {

	final static long serialVersionUID = (long) "VertexStack".hashCode();

	// array for all the z-slices; some Vertices null
	public Vertex[] vertices;
	
	public VertexStack(Vertex[] input) {
		vertices = input;
	}
	
	public Vertex[] vertices() {
		return vertices;
	}
	
	// the number of non-null Vertices in the stack
	public int size() {
		int size = 0;
		for (Vertex v : vertices) 
			if (v != null) size++;
		return size;
	}

	// returns all of the coordinates in a n x 2 array, where n = vertices.length
	public double[][] coordinateStack() {
		double[][] coordinates = new double[vertices.length][2];
		int i = 0;
		for (Vertex v : vertices) {
			if (v == null) {
				coordinates[i][0] = Double.NaN;
				coordinates[i][1] = Double.NaN;
			}
			else {
				coordinates[i][0] = v.coords()[0];
				coordinates[i][1] = v.coords()[1];
			}
			i++;
		}
		return coordinates;
	}
	
	public String toString() {
		return "VertexStack with " + size() + " vertices.";
	}
	
}
