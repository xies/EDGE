function handles = update_embryo(handles)
% if the user changes one of the base properties in the semiauto then we
% are in trouble because we don't want to lose all of the cellgraphs
% (obviously!!) so, here's what we do......
% we make a new embryo and copy over as many of the cellgraphs as possible.
% this means if you change to a smaller range of times, and then change
% back to the larger one, you've lost the CellGraphs from this process. but
% that's ok, the user has no reason to do this


% first, check that something is actually different........
% fields = changed_data_info(handles);
% if isempty(fields)  % if there are no changes
%     return
% end
% ************************!!!!!
%IN THE END WE WILL WANT THIS, BUT FOR NOW IT'S A GOOD HACK FOR FORCING IT
%TO TRRACK


readyproc(handles, 'tracking');
new_embryo = Embryo4D(handles.embryo, ...
    handles.info.start_time, handles.info.end_time, handles.info.master_time, ...
    handles.info.bottom_layer, handles.info.top_layer, handles.info.master_layer, ... 
    handles.info.tracking_area_change_Z, handles.info.tracking_layers_back_Z, ...
    handles.info.tracking_centroid_distance_Z / handles.info.microns_per_pixel, ...
    handles.info.tracking_area_change_T, handles.info.tracking_layers_back_T, ...
    handles.info.tracking_centroid_distance_T / handles.info.microns_per_pixel);
readyproc(handles, 'ready');

handles.embryo = new_embryo;
% save_embryo(handles);

% set the sliders to the reference image
handles = go_to_image(handles, handles.info.master_time, handles.info.master_layer);
