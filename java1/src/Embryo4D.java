import java.util.Vector;

public class Embryo4D implements java.io.Serializable {

	final static long serialVersionUID = (long) "Embryo4D".hashCode();
	
	public static final boolean DEBUG_MODE = true;
	
	// the maximum fraction by which the area of a Cell can change between slices
	private final double areaChangeMaxZ, areaChangeMaxT;
	// the number of layers to look back and still make a match
	private final int layersToLookBackZ, layersToLookBackT;
	// the maximum distance the centroids can differ between Cells to still make a match
	private final double centroidDistMaxZ, centroidDistMaxT;

	public int startTime, endTime, masterTime, bottomLayer, topLayer, masterLayer; // filename coords

	// all the data
	private CellGraph[][] cellGraphs;
	
	// the total number of (active) cells
	private int numCells;
	
	// has the Embryo been changed since the last save?
	public boolean changed;	
	
	public Embryo4D(int startTime,   int endTime,  int masterTime, 
					int bottomLayer, int topLayer, int masterLayer,
					double areaChangeMaxZ, int layersToLookBackZ, double centroidDistMaxZ, 
					double areaChangeMaxT, int layersToLookBackT, double centroidDistMaxT) {
		this.startTime = startTime;
		this.endTime = endTime;
		this.masterTime = masterTime;
		this.bottomLayer = bottomLayer;
		this.topLayer = topLayer;
		this.masterLayer = masterLayer;
		this.areaChangeMaxT = areaChangeMaxT;
		this.areaChangeMaxZ = areaChangeMaxZ;
		this.layersToLookBackT = layersToLookBackT;
		this.layersToLookBackZ = layersToLookBackZ;
		this.centroidDistMaxT = centroidDistMaxT;
		this.centroidDistMaxZ = centroidDistMaxZ;
		
		cellGraphs = new CellGraph[Math.abs(startTime-endTime)+1][Math.abs(bottomLayer-topLayer)+1];
		numCells = 0;
		changed = false;
		
		if (DEBUG_MODE && !isValid()) System.err.println("Error in Embryo4D:init!");
		assert(isValid());
	}
	public Embryo4D(Embryo4D oldEmbryo, int startTime,   int endTime,  int masterTime, 
			int bottomLayer, int topLayer, int masterLayer,
			double areaChangeMaxZ, int layersToLookBackZ, double centroidDistMaxZ, 
			double areaChangeMaxT, int layersToLookBackT, double centroidDistMaxT) {

		// call the normal constructor first
		this(startTime, endTime, masterTime, bottomLayer,  topLayer,  masterLayer,
				 areaChangeMaxZ,  layersToLookBackZ,  centroidDistMaxZ, 
				 areaChangeMaxT,  layersToLookBackT,  centroidDistMaxT);
		
		changed = true;
		
		// get whatever CellGraphs the old embryo has
		for (int i = 0; i < t(); i++) {
			for (int j = 0; j < z(); j++) {
				cellGraphs[i][j] = oldEmbryo.getCellGraph(unTranslateT(i), unTranslateZ(j));
				if (cellGraphs[i][j] != null)
					cellGraphs[i][j].setParent(this);
			}
		}
		
		// deactivate all the cells
		for (int i = 0; i < t(); i++)
			for (int j = 0; j < z(); j++)
				if (cellGraphs[i][j] != null)
					cellGraphs[i][j].deactivateAllCells();
			
		// perform tracking if it's full
		if (isFull()) 
			trackAllCells();
//			performTracking();
		
		if (DEBUG_MODE && !isValid()) System.err.println("Error in Embryo4D:init2!");
		assert(isValid());
	}
	
	// add a CellGraph a (t, z)
	// overwrite any existing CellGraph at that point
	public void addCellGraph(CellGraph cg, int t, int z) {
		int T = translateT(t);
		int Z = translateZ(z);
		cellGraphs[T][Z] = cg;
		cg.setParent(this);
		changed = true;
		
		if (isFull())
			trackAllCells();
//			performTracking();
	}
	
	// does this embryo contain all the necessary CellGraphs?
	private boolean isFull() {
		for (int i = 0; i < t(); i++)
			for (int j = 0; j < z(); j++)
				if (cellGraphs[i][j] == null)
					return false;
		return true;
	}
	// just a better name for accessing it from Matlab. but really
	// it does tracking whenever it's full, and therefore it's tracked whenever it's
	// full so they are equivalent. but it's weird to just call the function isTracked because
	// it is used as the criterion to determine whether to perform tracking. a statement like
	// if(isTracked()) performTracking();  would be weird.
	public boolean isTracked() {
		
		// actually, this function could really CHECK that the tracking works. this will be nCells times
		// faster than the tracking itself so it can be really fast- just a second. the reason is that the
		// cellAtPoint command takes at least order nCells time since it searches all cells. this way
		// we don't need to do that- we just make sure their centroids are inside of each other
		return isFull();
	}
	
	// deactivate all cells in a given direction. start at cell from but DO NOT INCLUDE
	// this cell.
	private void deactivateSingleCellZ(int c, int t, int from, int dir) {
		if (dir == 0) {
			if (cellGraphs[t][from].getCell(c) != null)
				cellGraphs[t][from].deactivateCell(c);
		}
		else if (dir == +1 || dir == -1) {
			for (int i = from + dir; i < z() && i >= 0; i += dir)
				if (cellGraphs[t][i].getCell(c) != null)
					cellGraphs[t][i].deactivateCell(c);
		}
	}
	private void deactivateSingleCellZ(int c, int t) {
		// at this rate could just loop through all z.......
		// in fact, all this from stuff is pointless
		int master_layer = translateZ(masterLayer);
//		deactivateSingleCellZ(c, t, master_layer,  0);
		deactivateSingleCellZ(c, t, master_layer, +1);
		deactivateSingleCellZ(c, t, master_layer, -1);
//		if (!isValid()) System.err.println("Error in Embryo4D:deactivateSingleCellZ!");
		assert(isValid());
	}
	// deactivate a given cell for all Z in time, starting at from and going in direction dir
	private void deactivateSingleCellTZ(int c, int from, int dir) {
		if (dir != 1 && dir != -1) return;
		for (int t = from + dir; t < t() && t >= 0; t += dir)
			for (int z = 0; z < z(); z++)
				if (cellGraphs[t][z].getCell(c) != null)
					cellGraphs[t][z].deactivateCell(c);
	}
	// does not deactivate at the master image-- i think this is the desired behavior(?)
	private void deactivateSingleCellTZ(int c) {
		int master_time = translateT(masterTime);
		int master_layer = translateZ(masterLayer);
		deactivateSingleCellZ(c, master_time, master_layer, +1);
		deactivateSingleCellZ(c, master_time, master_layer, -1);
		deactivateSingleCellTZ(c, master_time, +1);
		deactivateSingleCellTZ(c, master_time, -1);
//		if (!isValid()) System.err.println("Error in Embryo4D:decativateSingleCellTZ!");
		assert(isValid());
	}
	
	// spatially track the cell c at time t start at depth from in direction dir (+1 for up, -1 for down)
	private void trackSingleCellZ(int c, int t, int from, int dir) {
		if (dir != 1 && dir != -1) return;
			
		Cell cTrack = cellGraphs[t][from].getCell(c);
		// NOTE: this could be null for times other than the master time!!!
		// for these cells we simply don't do any tracking, unfortunately...?? I guess not.....
		// at least not for now. think about this..~~~
		if (cTrack == null)
			return;

		// track each layer (towards top)
		for (int i = from + dir; i < z() && i >= 0; i += dir) {				
			// for each Cell in the layer being tracked
			int j = i;
			// Tracking now looks FORWARD, not back!!
			for ( ; j < z() && j < i + layersToLookBackZ && j >= 0 && j > i - layersToLookBackZ; j += dir) {
				
				// the use of inactiveCellAtPoint means that you don't search for cell that have already been tracked
				// although actually that should probably never happen anyway.....
				// this is obsolete with the new system of more negative indices. ok... 
				// regular cellAtPoint would be fine too, it's the same thing
				Cell cMatchCandidate = cellGraphs[t][j].inactiveCellAtPoint(cTrack.centroid());
				if (cMatchCandidate == null) continue;
						
				/// this isCellMatch function contains the conditions for tracking
				if (isCellMatch(cTrack, cMatchCandidate, areaChangeMaxZ, centroidDistMaxZ)) {	
					cellGraphs[t][j].changeIndex(cMatchCandidate.index(), c);
					cTrack = cMatchCandidate; // now you track the next Cell to this one
					i = j;
					break;
				}

			} // look back some layers
			
			// if we reach here then we didn't find anything looking forward by that many layers.
			// in this case, we should give up.
			if (j == z() || j == i + layersToLookBackZ || j == -1 || j == i - layersToLookBackZ)
				break;
			
			
		} //
	}  // end of function

	private void trackSingleCellZ(int c, int t) {
		//here t is the array index form, it doesn't need to be "translated"
		
		// the actual array indices of the master values
//		int master_time = translateT(masterTime);
		int master_layer = translateZ(masterLayer);
		trackSingleCellZ(c, t, master_layer, +1);
		trackSingleCellZ(c, t, master_layer, -1);
//		if (!isValid()) System.err.println("Error in Embryo4D:trackSingleCellZ!");
//		assert(isValid());
	}
	
//	// (before doing this, need to activate everything in master layer)
//	// track a single cell from the master layer in both directions. as soon as you find a cell,
//	// you change it's index. for now this is just in Z
//	public void trackSingleCellZ(int c, int t) {	
//		// here t is the array index form, it doesn't need to be "translated"
//		
//		// the actual array indices of the master values
////		int master_time = translateT(masterTime);
//		int master_layer = translateZ(masterLayer);
//		
//		Cell cTrack = cellGraphs[t][master_layer].getCell(c);
//		// NOTE: this could be null for times other than the master time!!!
//		// for these cells we simply don't do any tracking, unfortunately...?? I guess not.....
//		// at least not for now. think about this..~~~
//		if (cTrack == null)
//			return;
//
//		// track each layer (towards top)
//		for (int i = master_layer + 1; i < z(); i++) {				
//			// for each Cell in the layer being tracked
//			int j = i;
//			for ( ; j < z() && j < i + layersToLookBackZ; j++) {
//				
//				// the use of inactiveCellAtPoint means that you don't search for cell that have already been tracked
//				// although actually that should probably never happen anyway.....
////				System.out.println(cellGraphs[t][j] + " " + t + " " + j);
////				System.out.println(cTrack);
//				Cell cMatchCandidate = cellGraphs[t][j].inactiveCellAtPoint(cTrack.centroidInt());
//				if (cMatchCandidate == null) continue;
//						
//				/// THIS isCellMatch function DEFINES THE CONDITIONS FOR TRACKING
//				if (isCellMatch(cTrack, cMatchCandidate, areaChangeMaxZ, centroidDistMaxZ)) {	
//					cellGraphs[t][j].changeIndex(cMatchCandidate.index(), c);
//					cTrack = cMatchCandidate; // now you track the next Cell to this one
//					i = j;
//					break;
//				}
//
//			} // look back some layers
//			
//			// if we reach here then we didn't find anything looking forward by that many layers.
//			// in this case, we should give up.
//			if (j == z() || j == i + layersToLookBackZ)
//				break;
//			
//		} // for each slice
//		
//		cTrack = cellGraphs[t][master_layer].getCell(c);
//
//		// Now track in the other direction (towards bottom)
//		for (int i = master_layer - 1; i >= 0; i--) {
//			int j = i;
//			for ( ; j >= 0 && j > i - layersToLookBackZ; j--) {	
//				
//				Cell cMatchCandidate = cellGraphs[t][j].inactiveCellAtPoint(cTrack.centroidInt());
//				if (cMatchCandidate == null) continue;
//				
//				/// THIS isCellMatch function DEFINES THE CONDITIONS FOR TRACKING
//				if (isCellMatch(cTrack, cMatchCandidate, areaChangeMaxZ, centroidDistMaxZ)) {	
//					cellGraphs[t][j].changeIndex(cMatchCandidate.index(), c);
//					cTrack = cMatchCandidate; // now you track the next Cell to this one
//					i = j;
//					break;
//				}
//											
//			} // look back some layers
//			
//			if (j == -1 || j == i - layersToLookBackZ)
//				break;
//			
//		} // for each slice
//
//		assert(isValid());
//	}
	
	// track cell c temporally start at time from in direction dir. tracking is always
	// carried out at the master layer
	private void trackSingleCellTZ(int c, int from, int dir) {
		if (dir != 1 && dir != -1) return;
		
		int master_layer = translateZ(masterLayer);	
		
		Cell cTrack = cellGraphs[from][master_layer].getCell(c);
		
		// who knows.......
		if (cTrack == null)
			return;
		
		// track each time point using mastser_time as a reference
		for (int i = from + dir; i < t() && i >= 0; i += dir) {  // i is a temporal variable
			int j = i;
			for ( ; j < t() && j < i + layersToLookBackT && j >= 0 && j > i - layersToLookBackT; j += dir) {
				Cell cMatchCandidate = cellGraphs[j][master_layer].inactiveCellAtPoint(cTrack.centroid());
				if (cMatchCandidate == null) continue;
				
				if (isCellMatch(cTrack, cMatchCandidate, areaChangeMaxT, centroidDistMaxT)) {
					cMatchCandidate.changeIndex(c);  // change the index to Z at the master layer 
					// (because cMatchCandidate is always from the master layer)
					trackSingleCellZ(c, j);  // this automatically changes the indices there
//					int cMatchCandidateIndex = cMatchCandidate.index();
//					for (int z = 0; z < z(); z++)
//						if (cellGraphs[j][z].getCell(cMatchCandidateIndex) != null)
//							cellGraphs[j][z].changeIndex(cMatchCandidateIndex, c);
					cTrack = cMatchCandidate;
					i = j;
					break;
				}
						
					
			} 
			
			// this code basically removes the looking back feature ... I don't think it's so good
			if (j == t() || j == i + layersToLookBackT || j == -1 || j == i - layersToLookBackT)
				break;
			
			
		} // for each time
		
	}
	
	private void trackSingleCellTZ(int c) {
		int master_time = translateT(masterTime);
		// track the cell spatially at the master time
		trackSingleCellZ(c, master_time);
		// track in time in both directions. this automatically tracks spatially at all other times
		trackSingleCellTZ(c, master_time, +1);
		trackSingleCellTZ(c, master_time, -1);
//		if (!isValid()) System.err.println("Error in Embryo4D:trackSingleCellTZ!");
//		assert(isValid());
	}

	
//	// assumes the cells are already spatially tracked!
//	public void trackSingleCellT(int c) {
//		// the actual array indices of the master values
//		int master_time = translateT(masterTime);
//		int master_layer = translateZ(masterLayer);
//				
//		Cell cTrack = cellGraphs[master_time][master_layer].getCell(c);
//		
//		// track each time point using mastser_time as a reference
//		for (int i = master_time + 1; i < t(); i++) {  // i is a temporal variable
//			int j = i;
//			for ( ; j < t() && j < i + layersToLookBackT; j++) {
//				Cell cMatchCandidate = cellGraphs[j][master_layer].inactiveCellAtPoint(cTrack.centroidInt());
//				if (cMatchCandidate == null) continue;
//				
//				if (isCellMatch(cTrack, cMatchCandidate, areaChangeMaxT, centroidDistMaxT)) {
//					int cMatchCandidateIndex = cMatchCandidate.index();
//					for (int z = 0; z < z(); z++)
//						if (cellGraphs[j][z].getCell(cMatchCandidateIndex) != null)
//							cellGraphs[j][z].changeIndex(cMatchCandidateIndex, c);
//					cTrack = cMatchCandidate;
//					i = j;
//					break;
//				}
//						
//					
//			} 
//			
//			if (j == t() || j == i + layersToLookBackT)
//				break;
//			
//		} // for each time
//		
//		cTrack = cellGraphs[master_time][master_layer].getCell(c);
//		
//		// round 2 -- track backward from master_time!!
//		for (int i = master_time - 1; i >= 0; i--) {  // i is a temporal variable
//			int j = i;
//			for ( ; j >= 0 && j > i - layersToLookBackT; j--) {
//				Cell cMatchCandidate = cellGraphs[j][master_layer].inactiveCellAtPoint(cTrack.centroidInt());
//				if (cMatchCandidate == null) continue;
//
//				if (isCellMatch(cTrack, cMatchCandidate, areaChangeMaxT, centroidDistMaxT)) {
//					int cMatchCandidateIndex = cMatchCandidate.index();
//					for (int z = 0; z < z(); z++)
//						if (cellGraphs[j][z].getCell(cMatchCandidateIndex) != null)
//							cellGraphs[j][z].changeIndex(cMatchCandidateIndex, c);
//					cTrack = cMatchCandidate;
//					i = j;
//					break;
//				}
//
//			} // look back some time steps
//			
//			if (j == -1 || j == i - layersToLookBackT)
//				break;
//			
//		} // for each time
//		
//		assert(isValid());
//	}

	public void trackAllCells() {
		int master_time = translateT(masterTime);
		int master_layer = translateZ(masterLayer);
				
		// first, deactivate all cells
		for (int t = 0; t < t(); t++)
			for (int z = 0; z < z(); z++)
				cellGraphs[t][z].deactivateAllCells();
		
		// deal with the master image (set the indices to all positive numbers)
		numCells = cellGraphs[master_time][master_layer].numCells();
		
		int i = 0;
		// for all cells in the master image
		for (Cell c : cellGraphs[master_time][master_layer].cells()) {
			// set the index of that cell to i, so that they get filled up from 1 to numCells
			cellGraphs[master_time][master_layer].changeIndex(c, ++i);
//			cellGraphs[master_time][master_layer].activateCell(c);
			// track
			trackSingleCellTZ(c.index()); // c.index() is just i
			
			// using this "activateCell" is a bit slower than using just the index i,
			// either way would be fine
		}
		
		
//		// spatially track for each time
//		for (int t = 0; t < t(); t++) {
//			int minIndex = CellGraph.minIndex(cellGraphs[t]);  // the minimum index of all cellgraphs (all z) at the given time
//			// that way all the cells in the master layer will have a smaller (thus unique) index
//			// now we don't use swapping anymore
//			Cell[] cells = cellGraphs[t][master_layer].cells();
//			for (Cell c : cells) {
//				if (t != master_time)
//					cellGraphs[t][master_layer].changeIndex(c, --minIndex);
//				// change the index to something even more negative
//				trackSingleCellZ(c.index(), t);
//			}
//		}
//
//		// temporally track
//		for (i = 1; i <= numCells; i++)
//			trackSingleCellT(i);
		if (DEBUG_MODE && !isValid()) System.err.println("Error in Embryo4D:trackAllCells!");
		assert(isValid());
	}

	
	

	// call this function if you modify the cell c at t, z (t and z already translated)
	// in necessary, retrack the cell
	private void retrackCell(int c, int t, int z) {
		if (!isTracked()) return;

		//		t = translateT(t);
//		z = translateZ(z);
		int master_time  = translateT(masterTime);
		int master_layer = translateZ(masterLayer);
				
		if (z == master_layer) {
			if (t == master_time) {  // z = master layer, t = master time
				// need to retrack the whole thing
				deactivateSingleCellTZ(c);
				trackSingleCellTZ(c);
			}
			else { // z = master layer,  t != master time
				deactivateSingleCellZ(c, t);
				deactivateSingleCellTZ(c, master_time, Misc.sign(t - master_time)); // for other times in that direction
				trackSingleCellZ(c, t);
				trackSingleCellTZ(c, master_time, Misc.sign(t - master_time));
			}
		}
		else {  // z != master layer
			// first, need to deactivate all cells in the direction you're tracking
			// in case they won't get tracked this time. note that this function does not deactivate
			// at the current depth z itself, it just prepares the others for tracking
			deactivateSingleCellZ(c, t, master_layer, Misc.sign(z - master_layer));
			// note this does not deactiviate the master_layer itself, which is good!
			
		
			// use signum so that if z > master_layer you track upwards, if z < master_layer you track downwards
			// for the 3rd argument, one might be tempted to save a bit more time and start from z instead
			// of master_layer (why re-track everything lower?). actually you'd need to start from 
			// something like z-layersToLookBackZ (or master_layer if you go past it when subtracked layers_to_look_backZ)
			// and I don't want to worry about that right now
			trackSingleCellZ(c, t, master_layer, Misc.sign(z - master_layer));
		}
		
		if (DEBUG_MODE && !isValid()) System.err.println("Error in Embryo4D:retrackCell!");
		assert(isValid());
	}
	
	
	// the point is to find the potential index of a new cell if it were to be tracked
	// to do this take a cell and "backtrack" to the master image.
	// if it runs into an active cell on the way, it immediately returns that index
	// otherwise, if it reaches the master image with no luck it returns 0.
	// ***actually we only need to go back one tracking step because everything
	// upstream is already tracked!!!!!!!!!
	private int backtrack(Cell cTrack, int t, int z) {
		Cell c = backtrackCell(cTrack, t, z);
		if (c == null) return 0;
		else return c.index();
	}
	
	
	/*
	 * the power of 1:1 tracking is that you can "run backwards down the ladder". for example, if you split a fat cell 
	 * into two (add edge), then you want to take both new cells, run down the Z-ladder to the master layer, then run 
	 * along the T-ladder to the master time, and now you have those indices because the indices are always guaranteed
	 *  to be there at the master image. then you just run back up (re-track) except you don't need to re-track 
	 *  everything, you just follow the rules in the re-track function as appropriate. the important part is that you 
	 *  got the cell indices necessary for doing that. I SEE. so,
if a cell is already active, you just use its index and go with it
if it's inactive, you run back to the master image and get its index. if you never reach the master image, you are 
hopeless and the re-track function returns (because you couldn't affect anything). if you are not hopeless then you
 just re-track as a said above.

now there's only one problem, and that's the layer_to_look_back thing. you can't exactly run backwards because you don't 
know where it jumped. if this value was exactly 2 it would be ok, but if its 3 or 2 or 1 then you get some ambiguities.
 perhaps just force it to be 2 (or 1?)????? 
wait, are there actually ambiguities? is there a case where you can get there in one direction but not along the other?
 well, actually, if tracking is 1:1, maybe not..É.

no, actually, it will be FINE with any layers to look back. that is an amazing result, rather surprising. see idea is
 that you don't always jump by layers_to_look_back, but rather you make the smallest jump possible <= that value. so 
 then you will make these same jumps in the other direction, aka you will always run right upto the edge where it fails 
 and then make the smallest jump. that's really cool.
 
 (23/3/10): some of the above is obsolete-- turns out you don't need 1:1 tracking to be able to run down the ladder, you 
 can just switch the order of the cells in the matching function. 
*/
	
	
	private Cell backtrackCell(Cell cTrack, int t, int z) {
		int master_layer = translateZ(masterLayer);
		int master_time  = translateT(masterTime);
		
		// first backtrack in z to the master layer
		if (z != master_layer) {
			// this is the OPPOSITE direction to that used in tracking. this is the "back" in backtrack..!!!
			int dir = Misc.sign(master_layer - z);
			
			for (int i = z + dir; i < z() && i < z + layersToLookBackZ && i >= 0 && i > z - layersToLookBackZ; i += dir) {
				// (back)tracking
				
				// in tracking we use inactiveCellAtPoint... here we use activeCellAtPoint
				// again, because we are going the right direction. if we find an inactive cell
				// that's OK, but it doesn't help us...
				Cell cMatchCandidate = cellGraphs[t][i].activeCellAtPoint(cTrack.centroid());
				if (cMatchCandidate != null) 				
					if (isCellMatch(cMatchCandidate, cTrack, areaChangeMaxZ, centroidDistMaxZ))
						// NOTE: cMatchCandidate and cTrack are switched!!!!!!!!!!!!!!!!!! this is crucial, so that
						// tracking actually doesn't need to be reversible! we just always compare them in a consistent
						// order, i.e., the one closest to the reference image is first, the farther one is second.
//						if (cMatchCandidate.isActive())
							return cMatchCandidate;
				
				if (i == master_layer)
					return null;
			}
		}
		else if (t != master_time) {  // at the master layer
			// moving towards master_time.....
			int dir = Misc.sign(master_time - t);
			for (int i = t + dir; i < t() && i < t + layersToLookBackT && i >= 0 && i > t - layersToLookBackT; i += dir) {
				
				// (back)tracking
				Cell cMatchCandidate = cellGraphs[i][master_layer].activeCellAtPoint(cTrack.centroid());
				if (cMatchCandidate != null) 				
					if (isCellMatch(cMatchCandidate, cTrack, areaChangeMaxT, centroidDistMaxT))
//						if (cMatchCandidate.isActive())   // implied by activeCellAtPoint
							return cMatchCandidate;
				
				if (i == master_time)
					return null;
			}
		}

		// if you never found anything, return null
		return null;
	}	
	

		
	
	// these are "high level" versions of add/remove/modify, where it also takes care
	// of any tracking issues. in other words, it calls retrackCell() in the appropriate
	// way depending on the nature of the cell
	public void addCell(Cell c, int t, int z) {
		t = translateT(t);
		z = translateZ(z);
		int master_time  = translateT(masterTime);
		int master_layer = translateZ(masterLayer);
		
		if (isTracked()) {
			if (t == master_time && z == master_layer) {
				cellGraphs[t][z].addCellActive(c);
				numCells++;
				retrackCell(c.index(), t, z);
			}
			else {
				cellGraphs[t][z].addCellInactive(c);
				int newInd = backtrack(c, t, z);
				if (newInd != 0) { // only want to add it and retrack if we found something to match it onto
					cellGraphs[t][z].changeIndex(c, newInd);
					retrackCell(newInd, t, z);
				}
			}
		}
		else {
			cellGraphs[t][z].addCellInactive(c);
		}
		changed = true;
		if (DEBUG_MODE && !isValid()) System.err.println("Error in Embryo4D:addCell!");
		assert(isValid());
	}
	public void removeCell(Cell c, int t, int z) {
		t = translateT(t);
		z = translateZ(z);
		int master_time  = translateT(masterTime);
		int master_layer = translateZ(masterLayer);
		
		if (CellGraph.isActive(c)) {
			if (t == master_time && z == master_layer) {
				deactivateSingleCellTZ(c.index());
				cellGraphs[t][z].removeCell(c);
				
				// want to keep things compact, so we move the last cell (index numCell)
				// to fill the spot that just opened by the deletion of c
				changeIndexAll(numCells, c.index());
				
				numCells--;
			}
			else {	
				int oldInd = c.index();
				cellGraphs[t][z].removeCell(c);
				retrackCell(oldInd, t, z);
			}
		}
		else {  // if inactive, just remove
			cellGraphs[t][z].removeCell(c);
		}
		changed = true;
		if (DEBUG_MODE && !isValid()) System.err.println("Error in Embryo4D:removeCell!");
		assert(isValid());
	}
	public void modifyCell(Cell c, int t, int z) {
		t = translateT(t);
		z = translateZ(z);
		
		if (isTracked()) {
			if (CellGraph.isActive(c))
				retrackCell(c.index(), t, z);
			else {  // if inactive
				int newInd = backtrack(c, t, z);
				if (newInd != 0) {
					cellGraphs[t][z].changeIndex(c, newInd);
					retrackCell(newInd, t, z);
				}
			}
		}
		// if the Embryo is not tracked, don't need to do anything
		changed = true;
		if (DEBUG_MODE && !isValid()) System.err.println("Error in Embryo4D:modifyCell!");
		assert(isValid());
	}

	
	// the number of time points
	public int t() {
		return cellGraphs.length;
	}
	// the number of space points
	public int z() {
		return cellGraphs[0].length;
	}
	// the number of Cells
	public int numCells() {
		return numCells;
	}
	public int Ys() {
		return cellGraphs[0][0].Ys;
	}
	public int Xs() {
		return cellGraphs[0][0].Xs;
	}

	// I want to keep the active Cell indices compact (from 1 to numCells)
	// so i need a changeIndex for all cells...
	private void changeIndexAll(int from, int to) {
		for (int t = 0; t < t(); t++)
			for (int z = 0; z < z(); z++)
				cellGraphs[t][z].changeIndex(from, to);
	}

	// the user gives time and space points on the interval [startTime endTime] 
	// and [bottomLayer topLayer] respectively. since the arrays in here go from
	// [0 t()-1] and [0 z()-1] respectively, we need to convert these coordinates.
	public int translateT(int t) {
		return Math.abs(t - startTime);
	}
	public int translateZ(int z) {
		return Math.abs(z - bottomLayer);
	}
	public int unTranslateT(int t) {
		return startTime + t * Misc.sign(endTime - startTime);
	}
	public int unTranslateZ(int z) {
		return bottomLayer + z * Misc.sign(topLayer - bottomLayer);
	}
	
	public CellGraph getCellGraph(int t, int z) {
		t = translateT(t);
		z = translateZ(z);
		if (t < 0 || t >= t() || z < 0 || z >= z()) return null;
		else 										return cellGraphs[t][z];
	}
	
	public Cell[] getCellStack(int c, int t) {
//		t = translateT(t);
//		Cell[] stack = new Cell[z()];
//		for (int i = 0; i < z(); i++)
//			stack[i] = cellGraphs[t][i].getCell(c);
//		return stack;
		return getCellStack(c, t, bottomLayer, topLayer + Misc.sign(topLayer-bottomLayer));
	}
	private Cell[] getCellStack(int c, int t, int zFrom, int zTo) {  // NOT including zTo
		t = translateT(t); 
		zFrom = translateZ(zFrom); zTo = translateZ(zTo);
		Cell[] stack = new Cell[Math.abs(zFrom-zTo)];
		int stackInd = 0;
		for (int i = zFrom; i != zTo; i+=Misc.sign(zTo-zFrom))
			stack[stackInd++] = cellGraphs[t][i].getCell(c);
		return stack;
	}
	public Cell[] getCellStackTemporal(int c, int z) {
//		z = translateZ(z);
//		Cell[] stack = new Cell[t()];
//		for (int i = 0; i < t(); i++)
//			stack[i] = cellGraphs[i][z].getCell(c);
//		return stack;
		return getCellStackTemporal(c, z, startTime, endTime + Misc.sign(endTime-startTime));
	}
	private Cell[] getCellStackTemporal(int c, int z, int tFrom, int tTo) {  // NOT including tTo
		z = translateZ(z); 
		tFrom = translateT(tFrom); tTo = translateT(tTo);
		Cell[] stack = new Cell[Math.abs(tFrom-tTo)];
		int stackInd = 0;
		for (int i = tFrom; i != tTo; i+=Misc.sign(tTo-tFrom))
			stack[stackInd++] = cellGraphs[i][z].getCell(c);
		return stack;
	}	
	
	public Cell getCell(int c, int t, int z) {
		// automatically translates the T, Z coordinate through getCellGraph fcn
		return getCellGraph(t, z).getCell(c);
	}
	public Cell[] getCells(int[] c, int t, int z) {
		Cell[] out = new Cell[c.length];
		for (int i = 0; i < c.length; i++)
			out[i] = getCell(c[i], t, z);
		return out;
	}
	
	// the highest tracked cell for a given cell index c at time t
	public double highestTracked(int c, int t) {
		if (!anyTracked(c, t)) return Double.NaN;
		for (int i = z()-1; i >= 0; i--)
			if (getCell(c, t, unTranslateZ(i)) != null)
				return unTranslateZ(i);
		return Double.NaN; // should never get here because of initial IF statement...
	}
	// the lowest tracked cell for a given cell index c at time t
	public double lowestTracked(int c, int t) {
		if (!anyTracked(c, t)) return Double.NaN;
		for (int i = 0; i < z(); i++)
			if (getCell(c, t, unTranslateZ(i)) != null)
				return unTranslateZ(i);
		return Double.NaN;
	}
	// was cell c tracked at all in this time point?
	public boolean anyTracked(int c, int t) {
		for (int i = 0; i < z(); i++)
			if (getCell(c, t, unTranslateZ(i)) != null)
				return true;
		return false;
	}
	public boolean isTracked(int c, int t, int z) {
		if (cellGraphs[t][z].getCell(c) != null) return true;
		else 									 return false;
	}
	
	
	
	// IMPORTANT. This function decides whether Cell c1 can be tracked to c2
	
	// ***************** the whole idea is that it is symmetric in the two cells !!!! ***
	// this is the 1:1 stuff that makes it all WORK (actually... don't need this property.. ha. 16/03/10)
	private boolean isCellMatch (Cell cTrack, Cell cMatchCandidate, 
			double area_change_max, double centroid_dist_max) {
		
		// STEP 0: Translate the Cell by the predicted amount, and then use the new translated cell
		// for all the below tests. 
		cMatchCandidate = predictLocation(cTrack, cMatchCandidate);
		
		
		// STEP 1: Centroid distance cannot be bigger than CENTROID_DIST_MAX
		if (Misc.distance(cTrack.centroid(), cMatchCandidate.centroid()) 
				> centroid_dist_max) return false;
		
		// STEP 2: Fractional area change cannot be bigger than AREA_CHANGE_MAX
//		double areaAverage = (cMatchCandidate.area() + cTrack.area()) / 2.0;
//		double areaRatio = Math.abs(cMatchCandidate.area() - cTrack.area()) / areaAverage;
//		if (areaRatio > area_change_max) return false;
		
		// STEP 2A: Overlap must be bigger than OVERLAP_MIN
		double min_overlap = area_change_max;
		if (overlapScore(cTrack, cMatchCandidate) < min_overlap) return false;
		// we use the area from cTrack because that is the trusted (i.e., already tracked) cell
		// and the area of cMatchCandidate could be something weird.
		// no, we use the max area. the problem with the above is that the new one can have a big portion
		// extra that doesnt exist in cTrack, and then nothing bad happens. to avoid this we use the max area,
		// since the max possible overlap area is that of the max area (i.e. when the two areas are equal and 
		// completely overlapping)
		
		// STEP 3: the two cells must contain each other's centroids
		if (!cMatchCandidate.containsPoint(cTrack.centroid())) return false;
		if (!cTrack.containsPoint(cMatchCandidate.centroid())) return false;
		// (note, the first of these is redundant, since that's how the match candidate is defined)
		// (but it's good to have them both here, since it ensures tracking is a 1:1 operation)
		// actually, there IS a reason to check both. if we are CHECKING the tracking (as 
		// opposed to doing it, then we want to check both of these conditions!!)
		// yeah.. so I dropped the 1:1 thing... but anyway... (16/03/10)
		
		return true;
	}
	private double overlapScore(Cell a, Cell b) {
		return a.overlapArea(b) / Math.max(a.area(), b.area());
	}
//	private double overlapScore2(Cell a, Cell b) {
//		return a.overlapArea(b) / Math.pow(Math.max(a.area(), b.area()), 2);
//	}
	// given a particular Cell cMatchCandidate, we want to translate it based on first order tilt (velocity) information
	// in order to let it line up with cTrack
	// this is accomplished by fitting a straight line through the centroids and using the slope to translate the
	// cell
	private Cell predictLocation(Cell cTrack, Cell cMatchCandidate) {
		final int MIN_PTS = 3;
		final int MAX_PTS = 7;
		
		// make a copy of the cell
		Cell c = new Cell(cMatchCandidate);
		
//		System.out.println(c.index() + "  t=" + c.t() + "  z=" + c.z());
		
		// get the centroids of all the cells from the master to this one
		double[][] centroids; double[] ztVals; double ztDelta;
		if (c.z() != masterLayer) {   // in this case we want a spatial stack
			centroids = Cell.centroidStack(getCellStack(cTrack.index(), cMatchCandidate.t(), masterLayer, c.z()));
			ztVals = new double[centroids.length];
//			for (int i = 0; i < translateZ(c.z()); i++) ztVals[i] = i;
			for (int i = 0; i < ztVals.length; i++) ztVals[i] = i;
			ztDelta = cMatchCandidate.z() - cTrack.z();
		}
		else {  // at masterLayer
			centroids = Cell.centroidStack(getCellStackTemporal(cTrack.index(), cMatchCandidate.z(), masterTime, c.t()));
			ztVals = new double[centroids.length];
//			System.out.println(ztVals.length + " " + c.t() + " " + translateT(c.t()));
			for (int i = 0; i < ztVals.length; i++) ztVals[i] = i;
			ztDelta = cMatchCandidate.t() - cTrack.t();
		}
		
		// find the number of non-NaN points, use only MAX_PTS of them starting from the end
		Vector<double[]> goodCentroids = new Vector<double[]>();
		Vector<Double  > goodZTvals    = new Vector<Double  >();
//		int ok = 0; 
//		for (int i = 0; i < centroids.length; i++)
//			if (centroids[i][0] != Double.NaN)
//				ok++;
		int ok = 0;
		for (int i = centroids.length-1; i >= 0; i--) {
			goodCentroids.add(centroids[i]);
			goodZTvals.add(ztVals[i]);
			if (centroids[i][0] != Double.NaN)
				ok++;
			if (ok >= MAX_PTS) 
				break;
		}
		// rewrite the centroids array so that in only includes at most MAX_PTS points
		double[][] centroidsTemp = new double[goodCentroids.size()][2];
		double[]   ztValsTemp    = new double[goodCentroids.size()];
		for (int i = 0; i < goodCentroids.size(); i++) {
			centroidsTemp[i] = goodCentroids.elementAt(i);
			ztValsTemp[i]    = goodZTvals.elementAt(i);
		}
		centroids = centroidsTemp;
		ztVals = ztValsTemp;
		
//		for (int i = 0; i < centroids.length; i++)
//		System.out.println("(" + centroids[i][0] + ", " + centroids[i][1] + ")   zt=" + ztVals[i] + "  i=" + i);
	
		
		// if the number of centroids is below MIN_PTS then don't do any translation
		if (ok >= MIN_PTS) {
			StraightLineFit fit = new StraightLineFit(centroids, ztVals);
			// in these units, each slice is separated by 1 unit. thus we multiply the slopes in 
			// x and y by 1 to get the translations in x and y.
			double[] translate = new double[2];
			translate[0] = -Math.abs(ztDelta) / fit.mY;
			translate[1] = -Math.abs(ztDelta) / fit.mX;
//			System.out.println(fit);
			// this is one approach. another approach is to use the whole fit to get a new location.
			// i am not sure which one i prefer. the translation method seems less precise but more robust (?)
			c.translate(translate);
			
//			System.out.println(cTrack);
//			System.out.println(cMatchCandidate);
//			System.out.println(c);
		}
		
		return c;
	}
	
	
	
	// this is some automatic error-correction code
	// that works as follow: for every inactive cell, check if it would be tracked to a cell right below it. 
	// (NOT POSSIBLE BY DEFINITION OF INACTIVE CELLS... (IT MEANS IMPOSE ONLY THE CENTROID INSIDE CONDITION,
	// NOT ALL THE OTHER TRACKING CONDITIONS)
	// if so, check if that cell has tracked to something on the current layer. if so, it means both cells should
	// have been tracked to the same one cell. thus they are probably the same cell, and we will merge
	// them using removeEdge(Cell a, Cell b) (assuming, of course, that they share an edge..!!)
	
	// ok...how do you decide where to add the edge? well, a cell (before relaxing edges) only has order 10 vertices,
	// so there are only << 100 possibilities to check. take the one which minimizes some metric. either centroid
	// distance or area change. let's say area change, because this centroid distance stuff might have to change
	// later with the "velocity" model if we decide to use that
	// no- area change is not good. we just need some "overlap" metric. this would be useful in several places,
	// such as for removing edges-- we have several options and want to consider the best one.
	
	public void autoRemoveEdges(int t, int z) {		
		CellGraph cg = getCellGraph(t, z);
		
		// for all inactive cells at this layer
		// do it this messy way because iterator has strange behavior when you mess around with things
		int[] inactiveCellIndices = cg.inactiveCellIndices();
		for (int loop = 0; loop < inactiveCellIndices.length; loop++) {
			int iThis = inactiveCellIndices[loop];
			// try merging with all of the neighbors (active and inactive, for now!)
			// for now only neighbors that share exactly 2 vertices (i.e., not more)
			Cell cThis = cg.getCell(iThis);
			
			
			// keep track of the potential merges and their "score"
			// then, at the end, pick the best one
			Vector<Cell[]> candidates = new Vector<Cell[]>();
			Vector<Double> candidatesOverlap = new Vector<Double>();
			
			
			for (int iNeigh : Cell.index(cg.cellNeighbors(iThis))) {
				
				Cell cNeigh = cg.getCell(iNeigh);
				
//				System.out.println(i);
//				System.out.println(cg.getCell(i));
//				System.out.println(neigh);
//				System.out.println(cg.getCell(neigh));
//				System.out.println();
				
//				// find the shared vertices
//				Vector<Vertex> verts = new Vector<Vertex>();
//				for (Vertex v : cg.vertices())
//					if (cg.getCell(i).containsVertex(v) && neighCell.containsVertex(v))
//						verts.add(v);
//				if (verts.size() != 2) continue;
//				// some adjacent cells might share more vertices if i refine edges
//				// this can be dealt with in theory, but for now it fails
//				// later should modify removeEdge to deal with this
				
				int newInd = cg.removeEdge(iThis, iNeigh);
//				System.out.println("newInd = " + newInd);
//				System.out.println(cg.getCell(newInd));
				
				// for now there are some problematic cases. for example, we cannot remove 
				// edges with >2 vertices 
				if (newInd == 0) continue;
				
				/* (16/03/10)
				 * this is not the right way to go about it, where i just pick the best neighbor to merge with  
				 * even in this case,it is not always an improvement!! after i find the best neighbor to merge with,
				 * i need to make sure that this merger actually _improves_ something. i.e., that the overlap with the 
				 * backtracked cell is actually MORE than it was before. NOTE:::: this only applies if i am an inactive
				 * cell merging with an _active_ cell. if i can merge with another inactive cell to make it tracked,
				 * then of course i should just go for it.
				 */ 
				 
				
				
				if (CellGraph.isActive(newInd)) {
//					System.out.println("Success! Created Cell " + newInd);
					
					Cell[] pair = new Cell[2];
					pair[0] = cThis;
					pair[1] = cNeigh;
					candidates.add(pair);
//					double overlap = cg.getCell(newInd).overlapArea(backtrackCell(cg.getCell(newInd), translateT(t), translateZ(z)));
//					overlap /= cg.getCell(newInd).area();
					candidatesOverlap.add(overlapScore(cg.getCell(newInd), backtrackCell(cg.getCell(newInd), translateT(t), translateZ(z))));
					
					// success
//					break; 
				}
//				else {
					
					// if I revert this with addEdge its too complicated. so i will manually
					// do this......
					// i want to remove and add the neigh in properly, because that might be active
					// but for "c", i don't care about adding it properly since it is inactive by definition
					// so i can just use the brute for adding mechanism. then i can also conveniently
					// set its index back to i. of course i could add the proper add and then change
					// index, i think either way would be fine
					
				removeCell(cg.getCell(newInd), t, z);
				cg.addCell(cThis, iThis);
				addCell(cNeigh, t, z);

				
			}  // for each neighbor
			
			// now pick the best cell (cell with highest overlap score)
			double max = 0;
			int maxInd = -1;
			for (int k = 0; k < candidates.size(); k++) {
				if (candidatesOverlap.elementAt(k) != Double.NaN && 
						candidatesOverlap.elementAt(k) > max) {
					max = candidatesOverlap.elementAt(k);
					maxInd = k;
				}
			}
			
			if (maxInd >= 0) {
				// (16/03/10) make sure that if the pair is active, then we actually want to merge~~~
				// by this i mean, either it is inactive, and that's fine, or it is active, and then we demand
				// that the "overlap score" is bigger now than it used to be.
				if (!candidates.elementAt(maxInd)[1].isActive() || 
						candidatesOverlap.elementAt(maxInd) > 
						overlapScore(candidates.elementAt(maxInd)[1], backtrackCell(candidates.elementAt(maxInd)[1], translateT(t), translateZ(z))))
				{
			
					// if you found anything, keep that Cell
			
					cg.removeEdge(candidates.elementAt(maxInd)[0], candidates.elementAt(maxInd)[1]);
				}
			}
			
		}   // for each Cell
		
		if (DEBUG_MODE && !isValid()) System.err.println("Error in Embryo4D: autoRemoveEdges!");
		assert(isValid());
	}
	
	public void autoAddEdges(int t, int z) {
		CellGraph cg = getCellGraph(t, z);
		
//		for (int H : cg.inactiveCellIndices())
//			System.out.println(H + ",  " + cg.getCell(H));
//		System.out.println("********");

		// for all inactive cells at this layer
		int[] inactiveCellIndices = cg.inactiveCellIndices();
		for (int loop = inactiveCellIndices.length - 1; loop >= 0; loop--) {
			int i = inactiveCellIndices[loop];			
			Cell c = cg.getCell(i);

// error might be happening because we take a copy of the array at the top, but when we re-add a cell
			// it might get a different index because inactives are not compact. that when we look
			// for the original it is nor there. is this right? why does it happen so rarely?
			// well, it's usually an EARLIEr one that gets lost...
//			for (int jj : cg.cellIndices())
//				if (cg.getCell(jj) == null)
//					System.err.println("Null cell in looper iteration " + i + " for Cell " + jj);
//			
//			for (int jj : inactiveCellIndices)
//				if (cg.getCell(jj) == null)
//					System.err.println("2 Null cell in looper iteration " + i + " for Cell " + jj);
			
			
//			System.out.println("****" + c);
					
			// keep track of the potential merges and their "score"
			// then, at the end, pick the best one
			Vector<Vertex[]> candidates = new Vector<Vertex[]>();
			Vector<Double> candidatesOverlap = new Vector<Double>();
				
			if (c == null) {
				System.err.println("WARNING (Embryo4D:autoAddEdges):null Cell during error correction! index = " + i);
				System.err.println("t = " + t + ", z = " + z);
			}
//			System.out.println(i);
			
			// for all unique pairs
			for (int vInd = 0; vInd < c.numV(); vInd++) {
				for (int wInd = vInd + 1; wInd < c.numV(); wInd++) {
					Vertex v = c.vertices()[vInd];
					Vertex w = c.vertices()[wInd];
					
					if (v == w) continue;
					if (cg.connected(v, w)) continue;
			
					// let this be handled by the correct cell-- don't let it happen for a cell
					// where the edge is not inside.....
					if (cg.cellAtPoint(Misc.midpoint(v.coords(), w.coords())) != c) continue;

					int[] newInds = cg.addEdge(v, w);
					// cases where the line between the new vertices is outside of the cell
					// and not inside any other cell (because that cell is at the edge of the embryo)
					
					if (newInds == null) continue;  
					// note in this case the edge was not added, so we don't need to worry
					
					int newIndex1 = newInds[0];
					int newIndex2 = newInds[1];
													
//					System.out.println("new inds = " + newIndex1 + ", " + newIndex2);
//					System.out.println(cg.getCell(newIndex1));
//					System.out.println(cg.getCell(newIndex2));
										
					if (CellGraph.isActive(newIndex1) && CellGraph.isActive(newIndex2)) {
						Vertex[] newPair = new Vertex[2];
						newPair[0] = v;
						newPair[1] = w;
						candidates.add(newPair);
						// get the overlap between the two new cells and the nearest backtracked cells
						// add them together for a total overlap score
						double overlap = overlapScore(cg.getCell(newIndex1), backtrackCell(cg.getCell(newIndex1), translateT(t), translateZ(z))) + 
										 overlapScore(cg.getCell(newIndex2), backtrackCell(cg.getCell(newIndex2), translateT(t), translateZ(z)));
							//cg.getCell(newIndex1).overlapArea(backtrackCell(cg.getCell(newIndex1), translateT(t), translateZ(z))) + 
								//		 cg.getCell(newIndex2).overlapArea(backtrackCell(cg.getCell(newIndex2), translateT(t), translateZ(z)));
						
						// the score will also be based on the smallest angle created by the new cells, because sometimes we 
						// get weird behavior that creates small angles
						double ANGLE_WEIGHTING = 1.0;  // note: angles are in radians
						overlap = overlap + ANGLE_WEIGHTING * (cg.getCell(newIndex1).vertexAngles()[0] + cg.getCell(newIndex2).vertexAngles()[0]); 
//						System.out.println(cg.getCell(newIndex1).vertexAngles()[0] + " and " + cg.getCell(newIndex1).vertexAngles()[1]);
						
						candidatesOverlap.add(overlap);
						
//						System.out.println("Success! Created cells.");
//						System.out.println("AAAAA " + cg.getCell(newIndex1).numV() + " " + cg.getCell(newIndex2).numV());
//						break;   // success
					}
										
					// remove that edge manually
					removeCell(cg.getCell(newIndex1), t, z);
					removeCell(cg.getCell(newIndex2), t, z);					
					cg.addCell(c, i);
					
					
				}  // each vertex
			}  // each vertex

			// now pick the best cell (cell with highest overlap score)
			double max = 0;
			int maxInd = -1;
			for (int k = 0; k < candidates.size(); k++) {
				if (candidatesOverlap.elementAt(k) != Double.NaN && 
						candidatesOverlap.elementAt(k) > max) {
					max = candidatesOverlap.elementAt(k);
					maxInd = k;
				}
			}
			// if you found anything, keep that Cell
			if (maxInd >= 0)
				cg.addEdge(candidates.elementAt(maxInd)[0], candidates.elementAt(maxInd)[1]);
			
			
		} // for each cell
		
		
		if (DEBUG_MODE && !isValid()) System.err.println("Error in Embryo4D: autoRemoveEdges!");
		assert(isValid());
		
	}
	
	
	private boolean isValid() {

		// this part is too slow...
//		for (int t = 0; t < t(); t++)
//			for (int z = 0; z < z(); z++)
//				if (cellGraphs[t][z] != null)
//					if (! cellGraphs[t][z].isValid())
//						return false;
		if (Math.abs(endTime - startTime) + 1 != t()) return false;
		if (Math.abs(bottomLayer - topLayer) + 1 != z()) return false;
		
		// it should be the case that if the cell is not tracked, then all cells are inactive
		if (! isTracked())
			for (int t = 0; t < t(); t++)
				for (int z = 0; z < z(); z++)
					if (cellGraphs[t][z] != null)
						if (cellGraphs[t][z].activeCells().length > 0) {
							System.err.println("Error: untracked Embryo contains active Cells");
							return false;
						}
		
		// the master image should contain only active cells
		int master_time = translateT(masterTime);
		int master_layer = translateZ(masterLayer);
		if (cellGraphs[master_time][master_layer] != null && isTracked()) { 
			if (cellGraphs[master_time][master_layer].inactiveCells().length > 0) {
				System.err.println("Error: reference CellGraph contains inactive Cells");
				return false;
			}
			
			// the set of active cell indices should be compact.
			int[] inds = cellGraphs[master_time][master_layer].activeCellIndices();
			for (int i = 0; i < numCells; i++) {
				if (inds[i] != i+1) {
					System.err.println("Error: indices not compact");
					return false;
				}
			}		
		}
		
		return true;
	}
	
	public String toString() {
		return "Embryo4D containing " + t() + " time points and " + z() + " z-slices";
	}
	

}
