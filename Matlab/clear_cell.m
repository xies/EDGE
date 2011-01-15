% Clear all graphics of a single cell, and removes it from the list of 
% selected cells. Only called by EDGE (not semiauto).

function handles = clear_cell(handles)

handles.activeCell = [];
handles.activeCellNeighbors = [];
set(handles.cell_text,'String', '-');

cla(handles.axes2);
set(handles.axes2, 'CameraPosition', 'default')
cla(handles.axes3);

% handles = slider_callbacks_draw_image_slice_dots_maingui(handles);