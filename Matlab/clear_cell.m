function handles = clear_cell(handles)
% clears a single cell. only called by EDGE (not semiauto);

handles.activeCell = [];
handles.activeCellNeighbors = [];
set(handles.cell_text,'String', '-');

cla(handles.axes2);
set(handles.axes2, 'CameraPosition', 'default')
cla(handles.axes3);

% handles = slider_callbacks_draw_image_slice_dots_maingui(handles);