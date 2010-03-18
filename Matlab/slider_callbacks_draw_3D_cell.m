function slider_callbacks_draw_3D_cell(handles, T, Z)



% this function is used both for making movies and for general use in EDGE.
% when making movies, we want to input an arbitrary T, Z, but in general we
% just use that from handles. therefore we can do it with 1 argument
% (normal) or 3 (movies)
if nargin == 1
    [T Z] = getTZ(handles);
    axes(handles.axes2);
    cla;  % clear axes
end

if isempty(handles.activeCell) || get(handles.button_manually_select_cells, 'Value')
    return;
end


DROPDOWN = get(handles.dropdown_measurements, 'Value');
drop_str = get(handles.dropdown_measurements, 'String');
MEASURE = drop_str{DROPDOWN};


 % draw 3d    
xlabel 'x (microns)'; ylabel 'y (microns)';
% title '3D reconstruction'
set(gca,'YDir','reverse')
hold on; 

% gets an array of cells making a "cell stack"
dx = handles.info.microns_per_pixel;
if get(handles.radiobutton_3d_spatial, 'Value')
    cell_stack = handles.embryo.getCellStack(handles.activeCell, T);
    zlabel 'z (microns)';
    dz = handles.info.microns_per_z_step;
    slice_highlight = abs(Z - handles.info.bottom_layer) + 1;
    % add 1 to match with matlab array indexing
else
    cell_stack = handles.embryo.getCellStackTemporal(handles.activeCell, Z);
    zlabel 't (minutes)';
    dz = handles.info.seconds_per_frame / 60;
    slice_highlight = abs(T - handles.info.start_time) + 1;
% add 1 to match with matlab array indexing
end



draw_cell_stack_highlight_z(cell_stack, T, slice_highlight, handles, dx, dz);


% draw all the neighbors as well
if ~isempty(handles.activeCellNeighbors) %&& get(handles.neighbors_3d, 'Value')
    for i = 1:length(handles.activeCellNeighbors)
        for j = 1:length(handles.activeCellNeighbors{i})
            cell_stack = handles.embryo.getCellStack(handles.activeCellNeighbors{i}(j), T);
            draw_cell_stack_highlight_z(cell_stack, T, Z,handles, dx, dz);
        end
    end
end


hold off;
axis equal;

% set(handles.text_readyproc, 'String', 'Ready');
% set(handles.text_readyproc, 'ForegroundColor', [0 1 0]);