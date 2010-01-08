import java.util.Arrays;
import java.util.Comparator;
import java.util.Vector;
import java.util.Set;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;

public class CellGraph implements java.io.Serializable {
	
	final static long serialVersionUID = (long) "CellGraph".hashCode();

	// the threshold for merging vertices, measured in pixels
	private static final double VERTEX_MERGE_DIST_THRESH_DEFAULT_VALUE = 3.0;
	
	// by default, don't do the min angle check
	private static final double VERTEX_MIN_ANGLE_DEFAULT_VALUE = 0.0;
	
	// the minimum number of neighbors a Cell and have
	private static final int MIN_CELL_NEIGHBORS = 2;  // do not allow cells or pairs of cells floating on their own
	
//	private Cell[] cells;  // no longer final because of mergeCells() aka removeEdge()
	private Map<Integer, Cell> cells = new TreeMap<Integer, Cell>();
	
	// size of the images
	public final int Ys;
	public final int Xs;
	
	// space and time coordinates of the CellGraph
	public final int t;
	public final int z;
	
	private Embryo4D parent;
	
	public static void help() {
		System.out.println("cells() - the cells");
		System.out.println("allVertices() - all the vertices in an (numVertices x 2) array");
		System.out.println("centroidCoords() - an (numCells x 2) array with the centroid coordinates");
		System.out.println("vertexCoords() - an (numVertices x 2) array with the vertex coordinates");
		System.out.println("numCells() - number of cells");
		System.out.println("numVertices() - number of vertices");
		System.out.println("cellsNeighboringVertex(Vertex) - all the cells neighboring the specified vertex");
		System.out.println("cellNeighbors(Cell) - returns the neighbors of that Cell in an array");
		System.out.println("connected(Cell, Cell) - are these two cells connected?");
		System.out.println("connected(Vertex, Vertex) - are these two vertices connected?");
		System.out.println("connectivityMatrixCell() - (numCells x numCells) cell connectivity matrix");
		System.out.println("connectivityMatrixVertex() - (numVertices x numVertices) vertex connectivity matrix");
		System.out.println("connectivityMatrixCellVertex() - an (numCells x numVertices) connectivity matrix");
		System.out.println("connectivityMatrixCellVertexOrdered() - a cell array of length numCells with the indices of the vertices of each cell");
		System.out.println("connectivityMatrixVertexInSameCell() - a connectivity matrix showing if vertices are in the same cell");
		System.out.println("draw() - returns a drawing of the CellGraph");
		System.out.println("drawCentroidConnectivity() - returns a drawing of the cell connectivity matrix");
		System.out.println("drawCell(int) - draws the individual cell the specified index");
		System.out.println("drawCentroids() - returns a drawing of all the centroids");
		System.out.println("drawVertices() - returns a drawing of all the vertices");
	}
	
		
	// ellipseProperties = YES; myosin = YES, distthresh = NO;
	public CellGraph(int[][] regions, int[][] centlist, int[][] vertlist, 
			int t, int z) {		
		this(regions, centlist, vertlist, t, z, VERTEX_MERGE_DIST_THRESH_DEFAULT_VALUE,
				VERTEX_MIN_ANGLE_DEFAULT_VALUE);
	}
	
	public CellGraph(int[][] regions, int[][] centlist, int[][] vertlist,
			int t, int z, double vertex_merge_dist_thresh, double vertex_min_angle) {		
		// assumes for regions: borders are -1, background is 0, cells are 1, 2, 3, ...
		// assumes for centlist/vertlist: two columns of coordinates, not two rows
		
		Ys = regions.length;
		Xs = regions[0].length;
		
		this.t = t;
		this.z = z;
		
		parent = null;
		
		// if it's an empty array
		if (centlist == null) {
//			cells = new Cell[0];
			return;
		}
		
		// temporary vertices (of type tempVertex) used until they are merged recursively
		Vector<TempVertex> tempVerts = new Vector<TempVertex>();
		
		// create the temporary vertices
		for (int i = 0; i < vertlist.length; i++) {
			Vertex v = new Vertex(vertlist[i][0], vertlist[i][1]);
			tempVerts.add(new TempVertex(v));
		}
		
		// remove those Vertices that border zero Cells
		
		// first, find which Cells each tempVertex touches
		for (int i = tempVerts.size() - 1; i >= 0; i--) {
			TempVertex v = tempVerts.elementAt(i);
				    
			// look at the 8 pixels around it
			for (int Y = (int) v.vertexCoords()[0] - 1; Y <= (int) v.vertexCoords()[0] + 1; Y++) {
				for (int X = (int) v.vertexCoords()[1] - 1; X <= (int) v.vertexCoords()[1] + 1; X++) { 
					// bounds checking
					if (Y <= 0 || Y > Ys || X <= 0 || X > Xs) continue;					
					
					// add the index of the Cell if the pixel is a Cell
					// v.cellNeighbors is a Set so duplicates will not be kept
					// EDITED SO THAT I ALSO KEEP ZEROS
//					if (regions[Y-1][X-1] > 0) v.addNeighbor(regions[Y-1][X-1] - 1);
					if (regions[Y-1][X-1] >= 0) v.addNeighbor(regions[Y-1][X-1] - 1);
	
				}
			}
			
			// keep it if there are > 0 neighbors
			if (v.neighbors().isEmpty()) tempVerts.remove(v);	
		}
		
		// this code goes with the fact that above I changed > to >=
		// basically sometimes I get a random small cell in the middle of a cell. this becomes
		// an extra vertex in the middle, and screws things up. this way if a vertex is in the middle
		// of a cell it only has one neighbor, whereas one on the edge that i actually want to keep
		// has 2 neighbors: that cell and the background (0). so i keep these first, and then remove them
		// (I want to get rid of vertices that are only touching one Cell AND are not at the edge)
		for (int i = tempVerts.size() - 1; i >= 0; i--) {
			TempVertex v = tempVerts.elementAt(i);

			if (v.numNeighbors() == 1) {
				tempVerts.remove(v);
			}
			// new: also, if a vertex has only two neighbors and is not along the edge, get rid of it
			// this prevents weird extra vertices from appearing in the middle
			else if (v.numNeighbors() == 2 && !v.containsNeighor(-1)) {
				tempVerts.remove(v);
			}
			else {
				// now can get rid of the background pixels
				v.removeNeighbor(-1);
			}
		}
	
		// calculate a distance matrix between the vertices
		// in preparation for merging the vertices together
		boolean[][] distMatrix = new boolean[tempVerts.size()][tempVerts.size()];
		
		for (int i = 0; i < tempVerts.size(); i++)
			for (int j = 0; j < tempVerts.size(); j++)
				distMatrix[i][j] = Misc.distance(tempVerts.elementAt(i).vertexCoords(), tempVerts.elementAt(j).vertexCoords()) <=
					vertex_merge_dist_thresh;
				
		// call the recursive merging function on each vertex
		Vector<TempVertex> mergedVerts = new Vector<TempVertex>();
		for (int i = 0; i < tempVerts.size(); i++) {
			// make sure the vertex is available
			if (!distMatrix[i][i]) continue;
			
			Vector<TempVertex> interVerts = new Vector<TempVertex>();
			TempVertex.vertexMergeRecursive(i, distMatrix, interVerts, tempVerts);
			
			mergedVerts.add(TempVertex.merge(interVerts));
		}
		
		Vector<Cell> initialCells = new Vector<Cell>();
		
		// create Cells
		for (int i = 0; i < centlist.length; i++) {
			Vector<Vertex> vertsOfCell = new Vector<Vertex>();
			
			// loop through all vertices
			for (int j = 0; j < mergedVerts.size(); j++) {
				// check each neighbor of that vertex
				for (int k : mergedVerts.elementAt(j).neighbors()) {		
					if (i == k) {
						vertsOfCell.add(mergedVerts.elementAt(j).vertex());
						break;
					}
				}	
			}
			
			// make sure that each Cell has at least 3 Vertices
			if (vertsOfCell.size() >= 3) {
				Vertex[] vertsOfCellArray = new Vertex[vertsOfCell.size()];
				vertsOfCell.toArray(vertsOfCellArray);
				
//				double myo;  // in case there's no myosin in the data set, just use myosinIntensity = 0
//				if (myosinIntensity == null)
//					myo = 0;
//				else
//					myo = myosinIntensity[i];
				Cell newCell = new Cell(centlist[i], vertsOfCellArray, this);
				initialCells.add(newCell);
			}
		}
		
		Cell[] initialCellsArray = new Cell[initialCells.size()];
		initialCells.toArray(initialCellsArray);
		
		Vector<Cell> finalCells = new Vector<Cell>();
	
		// now make sure each Cell has at least 2 Neighbors
		// and enforce the number of vertices if necessary
		for (Cell c : initialCellsArray)
			if (cellNeighbors(c, initialCellsArray).length >= MIN_CELL_NEIGHBORS)
				finalCells.add(c);
		
//		cells = new Cell[finalCells.size()];
//		finalCells.toArray(cells);
		
		for (int i = 0; i < finalCells.size(); i++)
			addCell(finalCells.get(i), -i - 1);  // strictly NEGATIVE INDICES --> INACTIVE

		/*
		// do the angle check  NOT SURE IF THIS IS RIGHT YET
		for (int i : cellIndices()) {
			Cell c = getCell(i);
			Vertex[] test = c.checkAngle(vertex_min_angle);
			if (test != null) {
//				System.out.println(test.length);
				removeCell(c);
			}
//			removeVertices(c.checkAngle(vertex_min_angle));
//			removeVertices(test);
			
//			if (c.area() == 0 || c.numV() < 3)
//				removeCell(c);
			
		}
		*/
			
		
		if (Embryo4D.DEBUG_MODE && !isValid()) System.err.println("Error in CellGraph:init!");
		assert(isValid());
	}
	
	public void setParent(Embryo4D parent) {
		this.parent = parent;
	}
	
	// returns the Cells
	public Cell[] cells() {
		Cell[] out = new Cell[cells.size()];
		cells.values().toArray(out);
		return out;
	}
	public Cell[] activeCells(Cell[] cells) {
		Vector<Cell> actives = new Vector<Cell>();
		for (Cell c : cells)
			if (isActive(c))
				actives.add(c);
		Cell[] out = new Cell[actives.size()];
		actives.toArray(out);
		return out;
	}
	public Cell[] activeCells() {
		return activeCells(cells());
	}
	public Cell[] inactiveCells(Cell[] cells) {
		Vector<Cell> actives = new Vector<Cell>();
		for (Cell c : cells)
			if (! isActive(c))
				actives.add(c);
		Cell[] out = new Cell[actives.size()];
		actives.toArray(out);
		return out;
	}
	public Cell[] inactiveCells() {
		return inactiveCells(cells());
	}
	public int[] cellIndices() {
		return Cell.index(cells());  // just the keySet
	}
	public int[] activeCellIndices() {
		return Cell.index(activeCells());
	}
	public int[] inactiveCellIndices() {
		return Cell.index(inactiveCells());
	}
	public boolean containsCell(int i) {
		if (getCell(i) != null) return true;
		else					return false;
	}
	
	
	// add Cell with index i
	public void addCell(Cell c, int i) {
		if (c == null) {
			System.err.println("Error in CellGraph:addCell-- adding null Cell!");
			return;
		}
		if (cells.containsValue(c)) 
			System.err.println("Warning: CellGraph/addCell(Cell, int) --- adding already present Cell " + c);
		if (cells.containsKey(i))
			System.err.println("Warning: CellGraph/addCell(Cell, int) --- overwriting cell " + i);		
		c.setIndex(i);
		cells.put(i, c);
//		if (Embryo4D.DEBUG_MODE && !isValid()) System.err.println("Error in CellGraph:addCell!");
		assert(isValid());
	}
	// add Cell and make its index the maximum available negative index
	public void addCellInactive(Cell c) {
		int i = 0;
		while (getCell(--i) != null);
		addCell(c, i);
	}
	// add Cell and make its index the minimum available positive index
	public void addCellActive(Cell c) {
		int i = 0;
		while (getCell(++i) != null);
		addCell(c, i);
	}
	public int minIndex() {
		int min = Integer.MAX_VALUE;
		for (int i : cells.keySet())
			min = Math.min(min, i);
		return min;
	}
	public static int minIndex(CellGraph[] cgs) {
		int min = cgs[0].minIndex();
		for (CellGraph cg : cgs)
			min = Math.min(min, cg.minIndex());
		return min;
	}
	public int maxIndex() {
		int max = Integer.MIN_VALUE;
		for (int i : cells.keySet())
			max = Math.max(max, i);
		return max;
	}
	public static int maxIndex(CellGraph[] cgs) {
		int max = cgs[0].maxIndex();
		for (CellGraph cg : cgs)
			max = Math.max(max, cg.maxIndex());
		return max;
	}
	
	// remove Cell with index i
	public Cell removeCell(int i) {
		return cells.remove(i);
	}
	public Cell removeCell(Cell c) {
		return removeCell(c.index());
	}
	public void activateCell(Cell c) {
		if (isActive(c)) return;
		addCellActive(removeCell(c));
	}
	public void activateCell(int i) {
		activateCell(getCell(i));
	}
	public void deactivateCell(Cell c) {
		if (!isActive(c)) return;
		addCellInactive(removeCell(c));
	}
	public void deactivateCell(int i) {
		deactivateCell(getCell(i));
	}
	public void deactivateAllCells() {
		Cell[] active = activeCells();
		for (Cell c : active)
			deactivateCell(c);
	}
	
	// get Cell with index i (returns null if this Cell does not exist)
	public Cell getCell(int i) {
		return cells.get(i);
	}
	// get Cells with indices i
	public Cell[] getCell(int[] inds) {
		if (inds == null) return null;
		Cell[] out = new Cell[inds.length];
		for (int i = 0; i < inds.length; i++)
			out[i] = getCell(inds[i]);
		return out;
	}
	
	// change the index of [the Cell with index i] from i to j
	public void changeIndex(int from, int to) {
		if (!containsCell(from)) return;
		
		// here I need to add a bit of a non-intuitive statement
		if (containsCell(to)) return;
		// in other words, you can't overwrite a cell with this function
		// the reasoning is as follows:
		/*
		when adding a new cell, you can get a situation where you overwrite an existing one
		 (in other words, this causes a problem because then tracking is not 1:1)
		here's an example:
		a cell in t= 3 tracks onto a cell (say, cell 5) in t=2, no problems
		then another cell also at t=3 (newly created by addEdge, for example) 
		does not track onto that cell at t=2, but it DOES track to cell 5 at t = 1, 
		so it gets the same index from backtrack and overwrites the (more rightful) cell 5 in t = 3. 
		so what's the deal here? you're always supposed to keep looking.
		we need to give priority to the DIRECT situation where they are at adjacent time points (i.e., the initial cell 5 at t=3)
		another reason why layers_to_look_back > 1 annoying...
		*/

		Cell removed = removeCell(from);
		addCell(removed, to);
	}
	public void changeIndex(Cell c, int to) {
		changeIndex(c.index(), to);
	}
//	public void changeIndexSwap(int from, int to) {
//		// here we do a special version in which
//		// we swap the indices if cell "to" already exists
//		// this is the if() statement. without the if statement it would
//		// just be regular changing of index where Cell "to" is overwritten it if exists
//		
//		Cell removed = cells.remove(from);  // could also call removeCell, doesn't matter
//		if (cells.containsKey(to))
//			changeIndex(to, from);
//		addCell(removed, to);
//	}
//	// change the index of [Cell c] to j
//	public void changeIndexSwap(Cell c, int to) {
//		changeIndexSwap(c.index(), to);
//	}
	// is the Cell c "active" ?
	public static boolean isActive(Cell c) {
		return isActive(c.index());
	}
	public static boolean isActive(int i) {
		if (i > 0) return true;
		else	   return false;
	}

	public int numCells() {
		return cells.size();
	}
	
	public int numActiveCells() {
		return activeCells().length;
	}
	public int numInactiveCells() {
		return inactiveCells().length;
	}
	
	public int numVertices() {
		return vertices().length;
	}
	public int numActiveVertices() {
		return activeVertices().length;
	}
	public int numInactiveVertices() {
		return inactiveVertices().length;
	}
	
	/// accessor methods ////
	private double[][] centroidCoords(Cell[] inputCells) {
		double[][] cents = new double[inputCells.length][2];
		for (int i = 0; i < inputCells.length; i++)
			cents[i] = inputCells[i].centroid();
		return cents;
	}
	// returns the coordinates of the centroids in an array
	public double[][] centroidCoords() {
		return centroidCoords(cells());
	}
	
	// returns the coordinates of the vertices in an array
	public double[][] vertexCoords() {
		return Vertex.coords(vertices());
	}
	
//	// returns the coordinates of the edges
//	public double[][][] activeEdges() {
//		return edges(activeCells());
//	}
//	public double[][][] inactiveEdges() {
//		return edges(inactiveCells());
//	}
//	public double[][][] edges() {
//		return edges(cells());
//	}
//	public double[][][] edges(Cell[] cells) {
//		// each each will go from a point (x1, y1) to a point (x2, y2)
//		Vector<Double> y1 = new Vector<Double>();
//		Vector<Double> x1 = new Vector<Double>();
//		Vector<Double> y2 = new Vector<Double>();
//		Vector<Double> x2 = new Vector<Double>();
//		
//		for (Cell c : cells) {
//			for (int i = 0; i < c.numV() - 1; i++) {
//				y1.add(c.vertices()[i].coords()[0]);
//				x1.add(c.vertices()[i].coords()[1]);
//				y2.add(c.vertices()[i+1].coords()[0]);
//				x2.add(c.vertices()[i+1].coords()[1]);
//			}
//			y1.add(c.vertices()[c.numV()-1].coords()[0]);
//			x1.add(c.vertices()[c.numV()-1].coords()[1]);
//			y2.add(c.vertices()[0].coords()[0]);
//			x2.add(c.vertices()[0].coords()[1]);
//		}
//		// the results will have 2 copies of each edge, but this is probably better (faster)
//		// than using a set. can change this easily;
//		
//		double[][][] out = new double[2][2][y1.size()];
//		double[] Y1 = new double[y1.size()];
//		double[] X1 = new double[y1.size()];
//		double[] Y2 = new double[y1.size()];
//		double[] X2 = new double[y1.size()];
//		for (int i = 0; i < y1.size(); i++) {
//			Y1[i] = y1.elementAt(i);
//			X1[i] = x1.elementAt(i);
//			Y2[i] = y2.elementAt(i);
//			X2[i] = x2.elementAt(i);
//		}
//		out[0][0] = X1;
//		out[0][1] = Y1;
//		out[1][0] = X2;
//		out[1][1] = Y2;
//		return out;
//	}
	
	
	/* the following functions MODIFY this CellGraph */
	
	public void removeCells(Cell[] toDelete) {
		for (Cell c : toDelete)
			parent.removeCell(c, t, z);
//			removeCell(c);
	}

	// these functions return the index of the new Cell if successful, or 0 if failure
	public int removeEdge(int c, int d) {
		return removeEdge(getCell(c), getCell(d));
	}
	// combine two cells into one by removing the edge between them
	// does nothing if the cells are not neighbors
	public int removeEdge(Cell c, Cell d) {
		if (!edgeConnected(c, d)) return 0;
		
		// find the two vertices that are shared by these cells
		// and then call removeEdge (Vertex, Vertex, Cell, Cell)
		Vector<Vertex> verts = new Vector<Vertex>();
		for (Vertex v : vertices())
			if (c.containsVertex(v) && d.containsVertex(v))
				verts.add(v);
		if (verts.size() != 2) return 0;
		
		Vertex v = verts.elementAt(0);
		Vertex w = verts.elementAt(1);
		return removeEdge(v, w, c, d);
	}
	public int removeEdge(Vertex v, Vertex w) {
		if (!connected(v, w)) return 0;
		
		// find the two cells that are associated with these two vertices,
		// and then call removeEdge (Vertex, Vertex, Cell, Cell)
		//  assumes there can only be two cells with this property
		Vector<Cell> cellVector = new Vector<Cell>();
		for (Cell can : cells())
			if (can.containsEdge(v, w))
				cellVector.add(can);
		if (cellVector.size() != 2) return 0;
		
		Cell c = cellVector.elementAt(0);
		Cell d = cellVector.elementAt(1);
		return removeEdge(v, w, c, d);
	}
	private int removeEdge(Vertex v, Vertex w, Cell c, Cell d) {
		// the vertices of the new Cell
		int newNumV = c.numV() + d.numV() - 2;  // -2 because they share v and w
		Vertex[] newVertArray = new Vertex[newNumV];

		// the following code is to preserve the counter-clockwise sorting
		// we don't want this to crash in the crazy case that one of the two cells is 
		// sorted clockwise by accident (otherwise we could just modify c.vertices(a, b) such that
		// it always started with the second one and went around back to the first). so i just assume
		// the first cell c is sorted properly (clockwise) and this way it will work when both are
		// ccw but when one is cw by accident it doesn't crash, it just might return a cw output (not that bad)
		int V = c.indexOf(v);
		int W = c.indexOf(w);
		if (W == V + 1 || (W == 0 && V == c.numV()-1)) {
			// need to swap them
			Vertex temp = v;
			v = w;
			w = temp;
		}
			
		Vertex[] cVerts = c.vertices(v, w); // from v, ending in w
		Vertex[] dVerts = d.vertices(w, v); // from w, ending in v
		
		// now that the vertices are sorted nicely, we just get all of them from c (except the last)
		// and then all of from d (except the last). we omit the last in each case so we don't count
		// the vertices v and w twice in the list for the new cell
		int i = 0;
		for (int j = 0; j < cVerts.length - 1; j++)
			newVertArray[i++] = cVerts[j];
		for (int j = 0; j < dVerts.length - 1; j++)
			newVertArray[i++] = dVerts[j];
		
		// use the special constructor that takes in a SORTED array of vertices
		Cell newCell = new Cell(newVertArray, this);

		
		// add in the new cell and remove the old ones
//		addCellInactive(newCell);
//		removeCell(c);
//		removeCell(d);
		
		// remove before adding for tracking reasons
		parent.removeCell(c, t, z);
		parent.removeCell(d, t, z);
		parent.addCell(newCell, t, z);
		return newCell.index();
	}
		
	// adds an edge between the vertices v and w. returns an array ret with ret[0] the index
	// of the cell removed, and ret[1] and ret[2] the indices of the two new cells
	public int[] addEdge(Vertex v, Vertex w) {
		// need to delete the cell in question, and add two cells
		
		// the first step is to find the cell that contains this edge
		// once we allow for splitting edges, there may be two cells that contain both of these 
		// (non-connected) vertices. thus, the strategy is to find the midpoint of the 
		// edge-to-be and check what cell it's inside. i guess this makes some assumptions (?)
		// about the shape of the cells but they are very reasonable ones
		Cell old = cellAtPoint(Misc.midpoint(v.coords(), w.coords()));
		if (old == null)
			return null;
		
		// the vertices of the new Cells
		Vector<Vertex> verts1 = new Vector<Vertex>();
		Vector<Vertex> verts2 = new Vector<Vertex>();
		boolean which = true;  // which cell are we currently adding vertices to
		
		// we go along toDelete's vertices and add them to the two groups
		// since this list is sorted, we keep assigning them to one of the new cells
		// until we hit v or w, which divides the new cells. then we add both
		// v and w to the new cells
		for (Vertex can : old.vertices()) {
			if (can == v || can == w) {
				which = !which;
				verts1.add(can);
				verts2.add(can);
			}
			else if (which)
				verts1.add(can);
			else
				verts2.add(can);	
		}
		// at the end of this algorithm, the vertex arrays come out SORTED for both
		// Cells. this is awesome, we use the special Cell constructor that inputs
		// a sorted vertex array and doesn't care about the centroid estimate.
		
		Vertex[] verts1Array = new Vertex[verts1.size()];
		Vertex[] verts2Array = new Vertex[verts2.size()];
		verts1.toArray(verts1Array);
		verts2.toArray(verts2Array);
		
		Cell new1 = new Cell(verts1Array, this);
		Cell new2 = new Cell(verts2Array, this);
		
//		removeCell(old);
//		addCellInactive(new1);
//		addCellInactive(new2);
		
		// need to remove first and then add (why? well, makes sense, don't want to overwrite in case of tracking)
		parent.removeCell(old, t, z);
		parent.addCell(new1, t, z);
		parent.addCell(new2, t, z);

			
		int[] ret = new int[2];
		ret[0] = new1.index();
		ret[1] = new2.index();
		return ret;
	}
	
	
//	public void activateCells(Cell[] inputCells) {
//		for (Cell c : inputCells)
//			activateCell(c);
//	}
	
	
	// make a cell out of the vertices V and add it
	public boolean addCell(Vertex[] verts) {
		// make sure there are >= 3 vertices
		if (verts.length < 3)
			return false;
		
		// make sure a cell with this set of vertices doesn't already exist
		// if all elements of A in B && and same length then we assume they are the same
		// because we assume no duplicate elements in either set
		HashSet<Vertex> vertsSet = new HashSet<Vertex>();
		for (Vertex v : verts)
			vertsSet.add(v);
		
		for (Cell c : cells()) {
			if (c.numV() == verts.length) {  // require that the number of Vertices are the same
				int i = 0;
				for ( ; i < c.numV(); i++) {
					if (!vertsSet.contains(c.vertices()[i]))
						break;  // break if that Vertex is not in the set
				}
				if (i == c.numV()) {// if they are all in the Set then this Cell already exists
					return false;
				}
			}
		}
		
		// make a crude estimate of the centroids by averaging the vertex positions
		int[] cent = new int[2];
		for (Vertex a : verts) {
			cent[0] = cent[0] + (int) a.coords()[0];
			cent[1] = cent[1] + (int) a.coords()[1];
		}
		cent[0] = cent[0] / verts.length;
		cent[1] = cent[1] / verts.length;

		// make the new cell
		Cell newCell = new Cell(cent, verts, this);

		// make sure it has positive area
		if (newCell.area() == 0)
			return false;
		
		parent.addCell(newCell, t, z);
//		addCellInactive(newCell);
		return true;
	}
	
	
	// supposed to be for the min angle check in the CG constuctor, but will
	// probably get rid of this...
	private void removeVerticesUntracked(Vertex[] verts) {
		if (verts == null) return;
		for (Vertex v : verts) {
			// for each Vertex we want to remove
			if (v == null) continue;
			
			for (Cell c : cellsNeighboringVertex(v)) {
				// for each Cell neighboring the Vertex in question
				if (c == null) continue;
				
				c.removeVertex(v);
				
				// if the cell becomes just a bunch of points of a line, remove it. 
				// i also use the redundant criterion numV < 3 because i'm worried about 
				// roundoff errors with getting an area of exactly zero, so now the worst
				// case scenario is just a cell consisting of >= 3 points on a line
				if (c.area() == 0 || c.numV() < 3) {
					removeCell(c);
					break;
				}
	
			}
		}
	}
	// remove the vertices verts
	// if it belongs to some Cells, need to fix that. 
	public void removeVertices(Vertex[] verts) {
		if (verts == null) return;
		for (Vertex v : verts) {
			// for each Vertex we want to remove
			if (v == null) continue;
			
			for (Cell c : cellsNeighboringVertex(v)) {
				// for each Cell neighboring the Vertex in question
				if (c == null) continue;
				
				c.removeVertex(v);
				
				// if the cell becomes just a bunch of points of a line, remove it. 
				// i also use the redundant criterion numV < 3 because i'm worried about 
				// roundoff errors with getting an area of exactly zero, so now the worst
				// case scenario is just a cell consisting of >= 3 points on a line
				if (c.area() == 0 || c.numV() < 3) {
					parent.removeCell(c, t, z);
					break;
//					removeCell(c);
				}
				else {
					parent.modifyCell(c, t, z);
				}
			
			}
		}
		if (Embryo4D.DEBUG_MODE && !isValid()) System.err.println("Error in CellGraph:removeVertices!");
		assert(isValid());
	}
	
	// split the edge between v and w to include another edge at coord
	public Vertex splitEdge(Vertex v, Vertex w, double[] coord) {
		if (!connected(v, w)) return null;
		
		Vertex newVert = new Vertex(coord[0], coord[1]);
		
		// find the two cells (at most 2) associated with this edge
		Cell c = null; 
		Cell d = null;
		for (Cell can : cells()) {
			if (can.containsVertex(v) && can.containsVertex(w)) {
				if (c == null)
					c = can;
				else
					d = can;
			}
		}
		
		c.splitEdge(v, w, newVert);
		parent.modifyCell(c, t, z);
		if (d != null) {
			d.splitEdge(v, w, newVert);
			parent.modifyCell(d, t, z);
		}
		
		if (Embryo4D.DEBUG_MODE && !isValid()) System.err.println("Error in CellGraph:splitEdge!");
		assert(isValid());
		return newVert;
	}
	
	public void moveVertex(Vertex v, double[] newcoords) {
		v.move(newcoords);
		for (Cell c : cellsNeighboringVertex(v))
			parent.modifyCell(c, t, z);
	}
	
	///////
	
	// refine the edges by splitting them
	public void refineEdges(boolean[][] bords, double maxAngle, double minAngle, double minEdgeLength) {
		Vertex[] vertices = vertices();
		for (Vertex v : vertices) {
			for (Vertex w : vertices) {
				if (connected(v, w)) {
					// if the edge length is *already* less than minEdgeLength, then by splitting it
					// things will only get worse. therefore we skip these cases
					if (Misc.distance(v.coords(), w.coords()) <= minEdgeLength) continue;
					// only do it once per edge
					refineEdge(v, w, bords, maxAngle, minAngle, minEdgeLength);
				}
			}
		}
		if (Embryo4D.DEBUG_MODE && !isValid()) System.err.println("Error in CellGraph:refineEdges!");
		assert(isValid());
	}
	// relax and edge by finding a pseudo-optimal vertex position for splitting the edge
	// we must create an angle on the interval [minAngle maxAngle]. the minimum prevents unrealistic
	// sharp corners, the maximum prevents the addition of vertices that don't capture any new curvature
	private boolean refineEdge(Vertex v, Vertex w, boolean[][] bords, double maxAngle, double minAngle, double minEdgeLength) {
		double[] cV = v.coords();
		double[] cW = w.coords();
		double[] midpt = Misc.midpoint(cV, cW);
		double slope = Misc.slope(cV, cW);
		double normal = -1/slope;
//		double edgeLength = Misc.distance(cV, cW);
		
		double[] r1 = new double[2];
		double[] r2 = new double[2];
		// doesn't have to be done this way, but I chose to move by 1 pixel each step
		double dx, dy;
		if (normal == 0) {
			dx = 1;
			dy = 0;
		}
		else if (slope == 0) {
			dx = 0;
			dy = 1;
		}
		else {
			dx =      1/Math.sqrt(normal*normal + 1);
		    dy = normal/Math.sqrt(normal*normal + 1);
		}
		
		for (int i = 0; i < Integer.MAX_VALUE; i++) {
			// move by 1 unit along the normal direction until less than 1 away from a border pixel
			r1[0] = midpt[0] + i * dy;
			r1[1] = midpt[1] + i * dx;
			
			r2[0] = midpt[0] - i * dy;
			r2[1] = midpt[1] - i * dx;
			
			double angle = Misc.angle(r1, cV, cW);  // angles should be the same for r1, r2

			// the termination statement of the loop
			if (angle < minAngle)
				return false;
//			if (angle > maxAngle)
//				continue;
			
			if (withinOne(r1, bords)) {
//			if (bords[(int)Math.round(r1[0])-1][(int)Math.round(r1[1])-1]) {
//				if ((Misc.distance(r1, cV) + Misc.distance(r1, cW) - edgeLength >= distThresh) &&
				 if (Misc.distance(r1, cV) >= minEdgeLength && Misc.distance(r1, cW) >= minEdgeLength &&
						 angle < maxAngle) {
					splitEdge(v, w, r1);
					return true;
				}
				else return false; // if one point doesn't pass the test, the other can't try
				 // for example, if nearest border has angle too big, then it shouldn't try
				 // the next very far away borders, the causes problems
			}
						
			if (withinOne(r2, bords)) {
//			if (bords[(int)Math.round(r2[0])-1][(int)Math.round(r2[1])-1]) {
//				if ((Misc.distance(r2, cV) + Misc.distance(r2, cW) - edgeLength >= distThresh) &&
				if (Misc.distance(r2, cV) >= minEdgeLength && Misc.distance(r2, cW) >= minEdgeLength && 
						angle < maxAngle) {
					splitEdge(v, w, r2);
					return true;
				}
				else return false;
			}
			
		}
		return false;
	}
	// is the point r at a TRUE in bords? checks a 1x2 square with r as the lower left corner
	private boolean withinOne(double[] r, boolean[][] bords) {
		int y = (int) Math.round(r[0]) -1;  // -1 because java uses 0-indexed arrays
		int x = (int) Math.round(r[1]) -1;
		
		if (y < 0 || x < 0 || y+1 >= bords.length || x+1 >= bords[0].length) return false;
		
//		if (bords[y][x] || bords[y+1][x] || bords[y][x+1] || bords[y+1][x+1]) return true;
		if (bords[y][x] || bords[y][x+1]) return true;
		
		else return false;	
	}	
	/*  end of modifying functions */
	
	
	public Vertex[] vertices(Cell[] cells) {
		Set<Vertex> vertSet = new HashSet<Vertex>();
		
		// loop through each cell and add all of its vertices to the Set
		for (Cell c : cells) 
			for (Vertex v : c.vertices())
				vertSet.add(v);
		
		// change the set into an array and return it
		Vertex[] verts = new Vertex[vertSet.size()];
		vertSet.toArray(verts);
		
		return verts;
	}
	public Vertex[] vertices() {
		return vertices(cells());
	}
	public Vertex[] activeVertices() {
		return vertices(activeCells());
	}
	public Vertex[] inactiveVertices() {
		return vertices(inactiveCells());
	}
	
	// ** more accessing methods **
	
	
	// returns the all the neighbors of c from 1st order to nth order
	// output array is in no particular order
	private Cell[] cellNeighborsFilled(Cell c, Cell[] cells, int n) {
		// if the cell itself is null (this shouldn't happen)
		if (c == null) return null;
		
		if (n == 0) {  // just return the input Cell c itself, packed in an array of size 1
			// this is the definition of 0th order-- we need this so that 1st order neighbors 
			// has this 0th order (aka the cell itself) removed
			Cell[] out = new Cell[1];
			out[0] = c;
			return out;
		}
		else if (n < 0) return null;  // should never happen!
		
		Set<Cell> cellSet = new HashSet<Cell>();
		cellSet.add(c);
		for (int i = 0; i < n; i++) {
			// for all the cells already in the set
			Cell[] tempSet = new Cell[cellSet.size()];
			cellSet.toArray(tempSet);
			for (Cell inSet : tempSet) {
				// find its neighbors and add them
				for (Cell all : cells) {
					if (connected(inSet, all))
						cellSet.add(all);
				}
			}
		} // this algorithm is a bit slow because for order i it computes the neighbors of ALL
		  // the lower orders every time (from 1 to i-2) instead of just i-1. could fix this at
		  // some point but i don't think it's a big deal
		
		// convert the Set into an array and return it
		Cell[] neighbors = new Cell[cellSet.size()];
		cellSet.toArray(neighbors);
		
		return neighbors;
	}
	// nth order neighbors (not filled in)
	public Cell[] cellNeighbors(Cell input, Cell[] cells, int n) {
		if (input == null) return new Cell[0];
		
		Set<Cell> cellSet = new HashSet<Cell>();
		for (Cell c : cellNeighborsFilled(input, cells, n))
			cellSet.add(c);
		for (Cell c : cellNeighborsFilled(input, cells, n-1))
			cellSet.remove(c);
		
		// convert the Set into an array and return it
		Cell[] neighbors = new Cell[cellSet.size()];
		cellSet.toArray(neighbors);
		
		// sorts the vertices in order of standard angle
        Comparator<Cell> byAngle = new ByAngle(input.centroid());
        //  sort starting from index 1 so that vertices[0] is first
		Arrays.sort(neighbors, 0, neighbors.length, byAngle);
		
		return neighbors;
	}
	public Cell[] cellNeighbors(Cell input, int n) {
		return cellNeighbors(input, cells(), n);
	}
	public Cell[] cellNeighbors(Cell input, Cell[] cells) {
		return cellNeighbors(input, cells, 1);
	}
	public Cell[] cellNeighbors(Cell input) {
		return cellNeighbors(input, cells());
	}
	public Cell[] cellNeighbors(int input) {
		return cellNeighbors(getCell(input));
	}
	public Cell[] cellNeighborsActive(Cell input, int n) {
		return cellNeighbors(input, activeCells(), n);
	}
	public Cell[] cellNeighborsActive(Cell input) {
		return cellNeighborsActive(input, 1);
	}
	
	// find the Cells touching the Vertex input
	public Cell[] cellsNeighboringVertex(Vertex input) {
		Vector<Cell> neighborVector = new Vector<Cell>();

		// for all cells, are they connected to input?
		for (Cell c : cells()) 
			if (c.containsVertex(input)) 
				neighborVector.add(c);

		// convert the Set into an array and return it
		Cell[] neighbors = new Cell[neighborVector.size()];
		neighborVector.toArray(neighbors);
		
		return neighbors;
	}

	// returns cellsNeighboringVertex(v).length for each Vertex
	public int[] numOfCellsNeighboringVertex() {
		Vertex[] input = vertices();
		int[] ret = new int[input.length];
		for (int i = 0; i < input.length; i++) 
			ret[i] = cellsNeighboringVertex(input[i]).length;
		return ret;
	}
	
	// are these two cells connected?
	// (takes constant time)
	public boolean connected(Cell a, Cell b) {
		for (Vertex v : a.vertices())
			for (Vertex w : b.vertices())
				if (v == w) return true;
		return false;
	}
	public boolean connected(int a, int b) {
		return connected(getCell(a), getCell(b));
	}
	
	// do these two cells share an edge?
	public boolean edgeConnected(Cell a, Cell b) {
		// these two cells must share exactly two vertices, and these vertices
		// must be connected
		Vertex[] shared = new Vertex[2];
		int share = 0;
		for (Vertex v : a.vertices()) {
			for (Vertex w : b.vertices()) {
				if (v == w) {
					if (share == 2) return false; // must not be more than two
					shared[share++] = v;
				}
			}
		}

		if (share == 2 && connected(shared[0], shared[1])) return true;
		else return false;
	}

	
//	// returns a connectivity matrix of the Cells
//	// (takes time proportional to #cells^2)
//	// note that the matrix will have some blank rows and columns
//	// because it is made to preserve the real indices of the cells
//	public boolean[][] connectivityMatrixCell() {
//		boolean[][] matrix = new boolean[numCells()][numCells()];
//		
//		for (int i = 0; i < numCells(); i++)
//			for (int j = 0; j < numCells(); j++)
//				if (connected(cells[i], cells[j])) matrix[i][j] = true;
//
//		return matrix;
//	}
//	
	
	public int[][] connectivityMatrixVertexSparse(Cell[] cells) {
		Vertex[] verts = vertices(cells);
		
		// we want to give each vertex a unique index for the duration of this call
		// therefore we use a Map with key = Vertex, value = index
		Map<Vertex, Integer> vertMap = new TreeMap<Vertex, Integer>();
		for (int i = 0; i < verts.length; i++)
			vertMap.put(verts[i], i+1);
		
		Vector<Integer> x = new Vector<Integer>();
		Vector<Integer> y = new Vector<Integer>();
	
		for (Cell c : cells) {
			for (int i = 0; i < c.numV() - 1; i++) {
				Vertex v1 = c.vertices()[i];
				Vertex v2 = c.vertices()[i+1];
				x.add(vertMap.get(v1));
				y.add(vertMap.get(v2));
			}
			Vertex v1 = c.vertices()[c.numV()-1];
			Vertex v2 = c.vertices()[0];
			x.add(vertMap.get(v1));
			y.add(vertMap.get(v2));
		}		
		
		int[][] matrix = new int[2][x.size()];
		for (int i = 0; i < x.size(); i++) {
			matrix[0][i] = x.elementAt(i);
			matrix[1][i] = y.elementAt(i);
		}

		return matrix;
	}
	public int[][] connectivityMatrixVertexSparse() {
		return connectivityMatrixVertexSparse(cells());
	}
	public int[][] connectivityMatrixVertexSparseActive() {
		return connectivityMatrixVertexSparse(activeCells());
	}
	public int[][] connectivityMatrixVertexSparseInactive() {
		return connectivityMatrixVertexSparse(inactiveCells());
	}
	
//	// returns a connectivity matrix of the Vertices
//	// (takes time proportional to #cells * #vertices^2)
//	// this should only take #vertices^2 time but that is impossible unless the Vertices
//	// know what their index is in the Vertices array. should consider this possibility!
//	// actually there is no Vertices array- just getVertices picks an arbitrary order
	public boolean[][] connectivityMatrixVertex(Cell[] inputCells) {
		Vertex[] verts = vertices();
		
		boolean[][] matrix = new boolean[verts.length][verts.length];
		for (int i = 0; i < verts.length; i++)
			for (int j = 0; j < verts.length; j++)
				if (connected(verts[i], verts[j], inputCells))
					matrix[i][j] = true;
		return matrix;
	}
	public boolean[][] connectivityMatrixVertex() {
		return connectivityMatrixVertex(cells());
	}

	
	// are these two vertices connected?
	// THEY MUST BE CONNECTED FOR ALL CELLS -- there are some weird cases (not sure how they arise)
	// where they are connected for one cell but not another
	// simple version would be:
	// for (Cell c : inputCells)
	// 	   if (c.connected(a, b)) return true;
	// return false
	public boolean connected(Vertex a, Vertex b, Cell[] inputCells) {
		int cellsContaining = 0;
		for (Cell c : inputCells) {
			if (c.containsVertex(a) && c.containsVertex(b)) {
				cellsContaining++;
				if (!c.connected(a, b)) return false;
			}
		}
		if (cellsContaining > 0) return true;
		else return false;
	}
	public boolean connected(Vertex a, Vertex b) {
		return connected(a, b, cells());
	}

	// are these two vertices both part of the same cell?
	private boolean inSameCell(Vertex a, Vertex b, Cell[] inputCells) {
		for (Cell c : inputCells)
			if (c.containsVertex(a) && c.containsVertex(b)) return true;
		return false;
	}
	public boolean inSameCell(Vertex a, Vertex b) {
		return inSameCell(a, b, cells());
	}
	
	////////
	// returns the Cell that the point coord lies in
	// returns null if it can't find one
	public Cell cellAtPoint(double[] coord) {
		return cellAtPoint(coord, cells());
	}
	public Cell activeCellAtPoint(double[] coord) {
		return cellAtPoint(coord, activeCells());
	}
	public Cell inactiveCellAtPoint(double[] coord) {
		return cellAtPoint(coord, inactiveCells());
	}
	private Cell cellAtPoint(double[] coord, Cell[] cells) {
		if (coord[0] <= 0 || coord[0] > Ys || coord[1] <= 0 || coord[1] > Xs)
			return null;
		for (Cell c : cells) 
			if (c.containsPoint(coord)) return c;
		return null;
	}
	
	
	// returns the Vertex nearest to the point coord and at most MIN_DIST_THRESH pixels away
	// returns null if it can't find one
	public Vertex vertexAtPoint(double[] coord, double MIN_DIST_THRESH) {
		if (coord[0] <= 0 || coord[0] > Ys || coord[1] <= 0 || coord[1] > Xs)
			return null;
		Vertex minDistV = null;
		double minDist = Double.POSITIVE_INFINITY;
		for (Vertex v : vertices()) {
			if (Misc.distance(v.coords(), coord) < minDist) {
				minDist = Misc.distance(v.coords(), coord);
				minDistV = v;
			}
		}
		// if the vertex is sufficiently close, allow it
		if (minDist < MIN_DIST_THRESH) return minDistV;
		else						   return null;
	}
	
	
	// draws the cells
	private double[][] draw(Cell[] cells) {
		double[][] image = new double[Ys][Xs];
		for (Cell c : cells) 
			c.draw(image);
		return image;
	}
	public double[][] draw() {
		return draw(cells());
	}	
	public double[][] drawActive() {
		return draw(activeCells());
	}
	public double[][] drawInactive() {
		return draw(inactiveCells());
	}
	
	
	// takes time Xs * Ys * nCells * C  where C is the point-in-polygon size and is a constant
	public double[][] drawRegions() {
		// first, draw the Cell regions. Then, draw the borders. That way the borders
		// can't have any holes in them
		double[][] image = new double[Ys][Xs];
		
		for (Cell c : cells()) {
			if (c == null) 
				continue;
				
			// draw in the Cell index at each point
			for (int y = 0; y < Ys; y++)
				for (int x = 0; x < Xs; x++)
					if (c.containsPoint(y, x))
						image[y][x] = c.index();
		}
		
		// draw the borders
		double borders[][] = draw();
		
		// set background to 0 and borders to -1, so that
		// cells can start from 1
		for (int y = 0; y < Ys; y++)
			for (int x = 0; x < Xs; x++)
				if (borders[y][x] > 0) 
				image[y][x] = -1;
		
		return image;
	}
	
	// draw all the Cell centroids as points
	public double[][] drawCentroids() {
		double[][] image = new double[Ys][Xs];
		for (Cell c : cells()) {
			int[] cent = c.centroidInt();
			image[cent[0] - 1][cent[1] - 1] = 1;
		}
		return image;
	}
	
	// draw all the Vertices as points
	public double[][] drawVertices() {
		double[][] image = new double[Ys][Xs];
		for (Vertex v : vertices()) {
			double[] coord = v.coords();
			image[(int) coord[0] - 1][(int) coord[1] - 1] = 1; 
		}
		return image;
	}
	
	public String toString() {
		return "CellGraph at t = " + t + " and z = " + z + 
		", of size " + Ys + " by " + Xs + " pixels (Ys by Xs) with " +
		numCells() + " Cells and " + numVertices() + " Vertices.";
	}
	
	public boolean isValid() {
		for (Cell c : cells()) {
			if (c == null || !c.isValid()) {
				System.err.println("Invalid Cell: " + c);
				return false;
			}
			for (Vertex v : c.vertices()) {
				double vY = v.coords()[0];
				double vX = v.coords()[1];
				if (vY <= 0 || vY > Ys || vX <= 0 || vX > Xs) {
					System.err.println("Invalid Vertex: " + v + " Coordinates out of bounds.");
					return false;
				}
			}
			
		}
		for (int i : cells.keySet()) {
			for (int j : cells.keySet()) {
				if (i != j && cells.get(i).equals(cells.get(j))) {
					System.err.println("Duplicate Cell: " + cells.get(i));
					return false;
				}
			}
		}			
		for (int i : cells.keySet()) {
			if (cells.get(i).index() != i) {
				System.err.println("Cell index does not match internal index for " + cells.get(i));
				return false;
			}
		}
		
		// make sure no cells have index exactly equal to zero
		for (int i : cells.keySet()) {
			if (i == 0) {
				System.err.println(cells.get(i) + " has index equal to 0 (not allowed).");
				return false;
			}
		}
	
		return true;
	}

	
/*************************************************************************/
	
	// for sorting by angle for the Vertices connected to a Cell
	private static class ByAngle implements Comparator<Cell> {
		private double[] origin;
        
        public ByAngle(double[] o) {
            origin = o;
        }
        
        // backwards so that a more negative angle is "greater"
        // this way the Vertices will be sorted clockwise
        public int compare(Cell a, Cell b) {
            double angleA = angle(origin, a);
            double angleB = angle(origin, b);
            if (angleA < angleB) return +1;
            if (angleA > angleB) return -1;
            else                 return  0;
        }
        
        // computes the angle that v makes with the point origin
        // returns the angle in the range [0, 2pi)
        // if v is at the origin 0.0 is returned
        private static double angle(double[] origin, Cell c) {
            double dx = c.centroid()[1] - origin[1];
            double dy = origin[0] - c.centroid()[0];
            double angle = Math.atan2(dy, dx);
            if (angle < 0) angle = 2*Math.PI + angle;
            return angle;
        }       
    }

}
