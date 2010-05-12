/**
 * A vertex in 2D space. Used to make up the {@link Cell} polygons.
 * 
 * @author Michael Gelbart
 *
 */
public class Vertex implements java.io.Serializable, Comparable<Vertex> {
	
	final static long serialVersionUID = (long) "Vertex".hashCode();

	// coordinates of the vertex in Matlab's order (y, x)
	private final double[] coords = new double[2];
	
	/** The standard constructor for class Vertex. */
	public Vertex(double y, double x) {
		coords[0] = y; 
		coords[1] = x;
	}
	
	/** Create a copy of input Vertex v.*/
	public Vertex(Vertex v) {
		coords[0] = v.coords[0];
		coords[1] = v.coords[1];
	}
	
	/** Get the coordinates of this Vertex. */
	public double[] coords() {
		return coords;
	}
	/** Get the coordinates of an array of Vertex objects. */
	public static double[][] coords(Vertex[] v) {
		if (v == null) return null;
		double[][] out = new double[v.length][2];
		for (int i = 0; i < v.length; i++)
			out[i] = v[i].coords;
		return out;
	}
	
	/** Create an array of new Vertex objects from a coordinate array. */
	public static Vertex[] createVertices(double[][] coords) {
		Vertex[] out = new Vertex[coords.length];
		for (int i = 0; i < coords.length; i++)
			out[i] = new Vertex(coords[i][0], coords[i][1]);
		return out;
	}
	
	/** Move the Vertex to a new location. */
	public void move(double[] newcoords) {
		if (newcoords.length != 2) return;
		coords[0] = newcoords[0];
		coords[1] = newcoords[1];
	} 
	/** Translate the position of the Vertex by delta. */
	public void translate(double[] delta) {
		if (delta.length != 2) return;
		coords[0] += delta[0];
		coords[1] += delta[1];
	}
	
	/** Draw a line connecting the two Vertices v and w. */
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
