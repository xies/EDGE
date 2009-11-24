function handles = maingui_initialize_neighbors(handles, val)
% initializes the neighbors for EDGE. this is called either when you change
% the neighbor number or when you click the mouse


handles.activeCellNeighbors = cell(val, 1);

if isempty(handles.activeCell)  % for example, if they turn off polygons
    return;
end

[T Z] = getTZ(handles);
for i = 1:val
    this_cell = handles.embryo.getCellGraph(T, Z).getCell(handles.activeCell);
    cellneighbors = ... 
        handles.embryo.getCellGraph(T, Z).cellNeighborsActive(this_cell, i);
    if isempty(cellneighbors)
        handles.activeCellNeighbors(i:end) = [];
        break;
    else
        handles.activeCellNeighbors{i} = Cell.index(cellneighbors);
    end
end