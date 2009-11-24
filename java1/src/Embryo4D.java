public class Embryo4D implements java.io.Serializable {

	final static long serialVersionUID = (long) "Embryo4D".hashCode();
	
	// the maximum fraction by which the area of a Cell can change between slices
	private final double areaChangeMaxZ, areaChangeMaxT;
	// the number of layers to look back and still make a match
	private final int layersToLookBackZ, layersToLookBackT;
	// the maximum distance the centroids can differ between Cells to still make a match
	private final double centroidDistMaxZ, centroidDistMaxT;

	public int startTime, endTime, masterTime, bottomLayer, topLayer, masterLayer;

	// all the data
	private CellGraph[][] cellGraphs;
	
	// the total number of (active) cells
	private int numCells;
	
	
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
		
		assert(isValid());
	}
	
	// add a CellGraph a (t, z)
	// overwrite any existing CellGraph at that point
	public void addCellGraph(CellGraph cg, int t, int z) {
		int T = translateT(t);
		int Z = translateZ(z);
		cellGraphs[T][Z] = cg;
		cg.setParent(this);
		
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
			for ( ; j < z() && j < i + layersToLookBackZ && j >= 0 && j > i - layersToLookBackZ; j += dir) {
				
				// the use of inactiveCellAtPoint means that you don't search for cell that have already been tracked
				// although actually that should probably never happen anyway.....
				// this is obsolete with the new system of more negative indices. ok... 
				// regular cellAtPoint would be fine too, it's the same thing
				Cell cMatchCandidate = cellGraphs[t][j].inactiveCellAtPoint(cTrack.centroidInt());
				if (cMatchCandidate == null) continue;
						
				/// THIS isCellMatch function DEFINES THE CONDITIONS FOR TRACKING
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
	}

	private void trackSingleCellZ(int c, int t) {
		//here t is the array index form, it doesn't need to be "translated"
		
		// the actual array indices of the master values
//		int master_time = translateT(masterTime);
		int master_layer = translateZ(masterLayer);
		trackSingleCellZ(c, t, master_layer, +1);
		trackSingleCellZ(c, t, master_layer, -1);
		assert(isValid());
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
				Cell cMatchCandidate = cellGraphs[j][master_layer].inactiveCellAtPoint(cTrack.centroidInt());
				if (cMatchCandidate == null) continue;
				
				if (isCellMatch(cTrack, cMatchCandidate, areaChangeMaxT, centroidDistMaxT)) {
					cMatchCandidate.changeIndex(c);  // change the index to Z at the master layer 
					// (because cMatchCandidate is always from the master layer)
					trackSingleCellZ(c, j);  // this automatically changes the indeces there
//					int cMatchCandidateIndex = cMatchCandidate.index();
//					for (int z = 0; z < z(); z++)
//						if (cellGraphs[j][z].getCell(cMatchCandidateIndex) != null)
//							cellGraphs[j][z].changeIndex(cMatchCandidateIndex, c);
					cTrack = cMatchCandidate;
					i = j;
					break;
				}
						
					
			} 
			
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
		assert(isValid());
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
				deactivateSingleCellTZ(c, master_time, (int) Math.signum(t - master_time)); // for other times in that direction
				trackSingleCellZ(c, t);
				trackSingleCellTZ(c, master_time, (int) Math.signum(t - master_time));
			}
		}
		else {  // z != master layer
			// first, need to deactivate all cells in the direction you're tracking
			// in case they won't get tracked this time. note that this function does not deactivate
			// at the current depth z itself, it just prepares the others for tracking
			deactivateSingleCellZ(c, t, master_layer, (int) Math.signum(z - master_layer));
			// note this does not deactiviate the master_layer itself, which is good!
			
		
			// use signum so that if z > master_layer you track upwards, if z < master_layer you track downwards
			// for the 3rd argument, one might be tempted to save a bit more time and start from z instead
			// of master_layer (why re-track everything lower?). actually you'd need to start from 
			// something like z-layersToLookBackZ (or master_layer if you go past it when subtracked LTLBZ)
			// and I don't want to worry about that right now
			trackSingleCellZ(c, t, master_layer, (int) Math.signum(z - master_layer));
		}
		
		
		
	}
	
	
	// the point is to find the potential index of a new cell if it were to be tracked
	// to do this take a cell and "backtrack" to the master image.
	// if it runs into an active cell on the way, it immediately returns that index
	// otherwise, if it reaches the master image with no luck it returns 0.
	// ***actually we only need to go back one tracking step because everything
	// upstream is already tracked!!!!!!!!!
	private int backtrack(Cell cTrack, int t, int z) {
		int master_layer = translateZ(masterLayer);
		int master_time  = translateT(masterTime);
		
		// first backtrack in z to the master layer
		if (z != master_layer) {
			// this is the OPPOSITE direction to that used in tracking. this is the "back" in backtrack..!!!
			int dir = (int) Math.signum(master_layer - z);
			
			for (int i = z + dir; i < z() && i < z + layersToLookBackZ && i >= 0 && i > z - layersToLookBackZ; i += dir) {
				// (back)tracking
				
				// in tracking we use inactiveCellAtPoint... here we use activeCellAtPoint
				// again, because we are going the right direction. if we find an inactive cell
				// that's OK, but it doesn't help us...
				Cell cMatchCandidate = cellGraphs[t][i].activeCellAtPoint(cTrack.centroidInt());
				if (cMatchCandidate != null) 				
					if (isCellMatch(cTrack, cMatchCandidate, areaChangeMaxZ, centroidDistMaxZ))
//						if (cMatchCandidate.isActive())
							return cMatchCandidate.index();
				
				if (i == master_layer)
					return 0;
			}
		}
		else if (t != master_time) {  // at the master layer
			// moving towards master_time.....
			int dir = (int) Math.signum(master_time - t);
			for (int i = t + dir; i < t() && i < t + layersToLookBackT && i >= 0 && i > t - layersToLookBackT; i += dir) {
				
				// (back)tracking
				Cell cMatchCandidate = cellGraphs[i][master_layer].activeCellAtPoint(cTrack.centroidInt());
				if (cMatchCandidate != null) 				
					if (isCellMatch(cTrack, cMatchCandidate, areaChangeMaxT, centroidDistMaxT))
//						if (cMatchCandidate.isActive())   // implied by activeCellAtPoint
							return cMatchCandidate.index();
				
				if (i == master_time)
					return 0;
			}
		}

		// if you never found anything, return 0
		return 0;
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
		assert(isValid());
	}
	public void removeCell(Cell c, int t, int z) {
		t = translateT(t);
		z = translateZ(z);
		int master_time  = translateT(masterTime);
		int master_layer = translateZ(masterLayer);
		
		if (cellGraphs[t][z].isActive(c)) {
			if (t == master_time && z == master_layer) {
				deactivateSingleCellTZ(c.index());
				cellGraphs[t][z].removeCell(c);
				numCells--;
				
				// want to keep things compact, so we move the last cell (index numCell)
				// to fill the spot that just opened by the deletion of c
				changeIndex(numCells, c.index());
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
			 
		assert(isValid());
	}
	public void modifyCell(Cell c, int t, int z) {
		t = translateT(t);
		z = translateZ(z);
		
		if (isTracked()) {
			if (cellGraphs[t][z].isActive(c))
				retrackCell(c.index(), t, z);
			else {  // if inactive
				int newInd = backtrack(c, t, z);
				if (newInd != 0) {
					cellGraphs[t][z].changeIndex(c, newInd);
					retrackCell(newInd, t, z);
				}
			}
		}
		// if the embryo is not tracked, don't need to do anything
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
	// so i need a changeindex for all cells...
	public void changeIndex(int from, int to) {
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
		return startTime + t * (int) Math.signum(endTime - startTime);
	}
	public int unTranslateZ(int z) {
		return bottomLayer + z * (int) Math.signum(topLayer - bottomLayer);
	}
	
	public CellGraph getCellGraph(int t, int z) {
		t = translateT(t);
		z = translateZ(z);
		if (t < 0 || t >= t() || z < 0 || z >= z()) return null;
		else 										return cellGraphs[t][z];
	}
	
	public Cell[] getCellStack(int c, int t) {
		t = translateT(t);
		Cell[] stack = new Cell[z()];
		for (int i = 0; i < z(); i++)
			stack[i] = cellGraphs[t][i].getCell(c);
		return stack;
	}
	public Cell[] getCellStackTemporal(int c, int z) {
		z = translateZ(z);
		Cell[] stack = new Cell[t()];
		for (int i = 0; i < t(); i++)
			stack[i] = cellGraphs[i][z].getCell(c);
		return stack;
	}
	
	public Cell getCell(int c, int t, int z) {
		// automatically translates the T, Z coordinate through getCellGraph fcn
		return getCellGraph(t, z).getCell(c);
	}
	
	// the highest tracked cell for a given cell index c at time t
	public double highestTracked(int c, int t) {
		if (!anyTracked(c, t)) return Double.NaN;
		int highest = -1;
		for (int i = z()-1; i >= 0; i--)
			if (getCell(c, t, unTranslateZ(i)) != null)
				highest = i;
		return unTranslateZ(highest);
	}
	// the lowest tracked cell for a given cell index c at time t
	public double lowestTracked(int c, int t) {
		if (!anyTracked(c, t)) return Double.NaN;
		int lowest = -1;
		for (int i = 0; i < z(); i++)
			if (getCell(c, t, unTranslateZ(i)) != null)
				lowest = i;	
		return unTranslateZ(lowest);
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
	// this is the 1:1 stuff that makes it all WORK
	private boolean isCellMatch (Cell cTrack, Cell cMatchCandidate, 
			double area_change_max, double centroid_dist_max) {
		// STEP 1: Centroid distance cannot be bigger than CENTROID_DIST_MAX
		if (Misc.distance(cTrack.centroid(), cMatchCandidate.centroid()) 
				> centroid_dist_max) return false;
		
		// STEP 2: Fractional area change cannot be bigger than AREA_CHANGE_MAX
		double areaAverage = (cMatchCandidate.area() + cTrack.area()) / 2.0;
		double areaRatio = Math.abs(cMatchCandidate.area() - cTrack.area()) / areaAverage;
		if (areaRatio > area_change_max) return false;
		
		// STEP 3: the two cells must contain each other's centroids
		if (!cMatchCandidate.containsPoint(cTrack.centroidInt())) return false;
		if (!cTrack.containsPoint(cMatchCandidate.centroidInt())) return false;
		// (note, the first of these is redundant, since that's how the match candidate is defined)
		// (but it's good to have them both here, since it ensures tracking is a 1:1 operation)
		// actually, there IS a reason to check both. if we are CHECKING the tracking (as 
		// opposed to doing it, then we want to check both of these conditions!!)
		
		return true;
	}
	
	
	
	// this is some automatic error-correction code
	// that works as follow: for every inactive cell, check if it would be tracked to a cell right below it. 
	// if so, check if that cell has tracked to something on the current layer. if so, it means both cells should
	// have been tracked to the same one cell. thus they are probably the same cell, and we will merge
	// them using removeEdge(Cell a, Cellb) (assuming, of course, that they share an edge..!!)
	
	// this could have been done during the tracking stage, actually. check all cells that track to a given cell,
	// merge them
	
	// only allow this if the new cell follows the tracking criteria. 
	
	// complication: if you find a pair of two cells that are both not tracked AND neighboring, consider if their
	// combination could be tracked to something. GOOD! 
	
	
	// 00
	// okay, that is like automatic "add edge". what about automatic remove edge? same thing! if you have a cell
	// that is untracked, then see if you can add an edge and make it work. --> this is all because of area thresholds, etc
	
	// ok...how do you decide where to add the edge? well, a cell (before relaxing edges) only has order 10 vertices,
	// so there are only << 100 possibilities to check. take the one which minimizes some metric. either centroid
	// distance or area change. let's say area change, because this centroid distance stuff might have to change
	// later with the "velocity" model if we decide to use that
	
	// do we want to do this in space or time? or is that the point? should we do it spatially after
	// spatial tracking, and then temporally after temporal tracking? then it would be done during tracking
	
	// IDEA: maybe we should ONLY allow these things to act on inactive (or 1 active and 1 inactive in the first case)
	// that way it preserves tracking and we don't need to re-track. like even for the manual changes... well.......
	// not sure about this..~~
	
	public void autoAddEdges() {
		for (int i = translateZ(masterLayer) + 1; i < z(); i++) {		
			
			
		}
		for (int i = translateZ(masterLayer) - 1; i > 0; i--) {		
			
			
		}
	}
	
	
	
	private boolean isValid() {
		for (int t = 0; t < t(); t++)
			for (int z = 0; z < z(); z++)
				if (cellGraphs[t][z] != null)
					if (! cellGraphs[t][z].isValid())
						return false;
		if (Math.abs(endTime - startTime) + 1 != t()) return false;
		if (Math.abs(bottomLayer - topLayer) + 1 != z()) return false;
		
		// it should be the case that if the cell is not tracked, then all cells are inactive
		if (! isTracked())
			for (int t = 0; t < t(); t++)
				for (int z = 0; z < z(); z++)
					if (cellGraphs[t][z].activeCells().length > 0)
						return false;
		
		// the master image should contain only active cells
		int master_time = translateT(masterTime);
		int master_layer = translateZ(masterLayer);
		if (cellGraphs[master_time][master_layer].inactiveCells().length > 0)
			return false;
		
		// the set of active cell indices should be compact.
		int[] inds = cellGraphs[master_time][master_layer].activeCellIndices();
		for (int i = 0; i < numCells; i++) {
			if (inds[i] != i+1) {
				System.err.println("Error: indices not compact");
				return false;
			}
		}
		
		return true;
	}
	
	public String toString() {
		return "Embryo4D containing " + t() + " time points and " + z() + " slices";
	}
	
	
/*************************************************************************/
//	
//	// for sorting by angle for the Vertices connected to a Cell
//	private static class ByAngle implements Comparator<Cell4D> {
//		private double[] origin;
//		private int t;
//		private int z;
//        
//        public ByAngle(double[] o, int t, int z) {
//            origin = o;
//            this.t = t;
//            this.z = z;
//        }
//        
//        // backwards so that a more negative angle is "greater"
//        // this way the Vertices will be sorted clockwise
//        public int compare(Cell4D a, Cell4D b) {
//            double angleA = angle(origin, a.cells()[t].cells()[z]);
//            double angleB = angle(origin, b.cells()[t].cells()[z]);
//            if (angleA < angleB) return +1;
//            if (angleA > angleB) return -1;
//            else                 return  0;
//        }
//        
//        // computes the angle that v makes with the point origin
//        // returns the angle in the range [0, 2pi)
//        // if v is at the origin 0.0 is returned
//        private static double angle(double[] origin, Cell c) {
//        	return Misc.angle(origin, c.centroid());
//        }       
//    }


}
