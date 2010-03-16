public class Vertex implements java.io.Serializable, Comparable<Vertex> {
	
	final static long serialVersionUID = (long) "Vertex".hashCode();

	// coordinates of the vertex in Matlab's order (y, x)
	private final double[] coords = new double[2];
	
	public static void help() {
		System.out.println("coords() - the coordinates of the Vertex");
		System.out.println("move(int[] newcoords) - moves the Vertex to newcoords");
	}
	
	public Vertex(double y, double x) {
		coords[0] = y; 
		coords[1] = x;
	}
	
	public Vertex(Vertex v) {
		coords[0] = v.coords[0];
		coords[1] = v.coords[1];
	}
	
	// return the coordinates in an array
	public double[] coords() {
		return coords;
	}
	public static double[][] coords(Vertex[] v) {
		if (v == null) return null;
		double[][] out = new double[v.length][2];
		for (int i = 0; i < v.length; i++)
			out[i] = v[i].coords;
		return out;
	}
	
	// create a bunch of Vertex objects from the coordinates
	public static Vertex[] createVertices(double[][] coords) {
		Vertex[] out = new Vertex[coords.length];
		for (int i = 0; i < coords.length; i++)
			out[i] = new Vertex(coords[i][0], coords[i][1]);
		return out;
	}
	
	public void move(double[] newcoords) {
		coords[0] = newcoords[0];
		coords[1] = newcoords[1];
	}
	
	// draw a line connecting the two vertices v and w
	public static void drawConnection(double[][] image, Vertex v, Vertex w) {
		Misc.drawLine(image, v.coords, w.coords);
	}
	
	public String toString() {
		return "Vertex at (" + coords[0] + ", " + coords[1] + ")";
	}

	public int compareTo(Vertex that) {
		if (this.coords[0] > that.coords[0]) return +1;
		if (this.coords[0] < that.coords[0]) return -1; 
		if (this.coords[1] > that.coords[1]) return +1;
		if (this.coords[1] < that.coords[1]) return -1;
		else 								 return 0;
	}
	
}
