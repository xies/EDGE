% Used for changing data sets in semiauto. Only does tasks specific to 
% semiauto, i.e., those that don't apply in EDGE such as tracking parameters.

function handles = clear_data_set_semiauto(handles, data_set)
 
readyproc(handles, 'loading')

handles.is_semiauto = 1;
handles.is_maingui = 0;

%%%% get data info
info = read_data_info(data_set);


if info.notfound  % if this is a new data set, never been imported
    msgbox(strcat('Info for data set "', data_set, ...
        '" does not exist. There has been an error in the import.', ...
        'Check DATA_INFO.csv to see that an entry for this data set was created correctly.'), ...
        'Loading failed', 'error'); 
    return;
end

handles.data_set = data_set;
handles.info = info;

% set the source paths
handles = handles_set_program_dir(handles);
handles = handles_set_src_paths(handles);
handles = handles_set_tempsrc_paths(handles);

handles = clear_data_set_both(handles); % takes care of most things

handles = initialize_tz_sliders(handles);


% load the embryo object
if ~exist(fullfile(handles.tempsrc.parent, 'embryo_data.mat'), 'file')
%         handles.embryo = Embryo4D(...
%         handles.info.start_time, handles.info.end_time, handles.info.master_time, ...
%         handles.info.bottom_layer, handles.info.top_layer, handles.info.master_layer, ... 
%         handles.info.tracking_area_change_Z, handles.info.tracking_layers_back_Z, ...
%         handles.info.tracking_centroid_distance_Z / handles.info.microns_per_pixel, ...
%         handles.info.tracking_area_change_T, handles.info.tracking_layers_back_T, ...
%         handles.info.tracking_centroid_distance_T / handles.info.microns_per_pixel);
%         save_embryo(handles);
%     
    msgbox('Cannot find DATA_SEMIAUTO embryo.mat file.', 'Failure', 'error');
    return
end
handles.embryo = [];
load(fullfile(handles.tempsrc.parent, 'embryo_data'));
handles.embryo = embryo4d;

% take out the temporal tracking parameters for fixed data sets

if handles.fixed
    set(handles.text_tracking_parameters_T, 'Visible', 'off');
    set(handles.info_text_tracking_area_change_T, 'Visible', 'off');
    set(handles.info_text_tracking_layers_back_T, 'Visible', 'off');
    set(handles.info_text_tracking_centroid_distance_T, 'Visible', 'off');
else
    set(handles.text_tracking_parameters_T, 'Visible', 'on');
    set(handles.info_text_tracking_area_change_T, 'Visible', 'on');
    set(handles.info_text_tracking_layers_back_T, 'Visible', 'on');
    set(handles.info_text_tracking_centroid_distance_T, 'Visible', 'on');    
end


if exist(fullfile(handles.src.parent, 'embryo_data.mat'), 'file')
    exported_text_controller(handles, 'exported');
else
    exported_text_controller(handles, 'not exported');
end
set(handles.cbox_inactive, 'Value', 1);

% clear the current cells
handles.activeCell = [];
handles.activeVertex = [];
handles.activeAdjustment = [];

% this is for automatically error correct *some* images-- small detail
handles.some_auto_range = [];

% initialize the vec buttons
handles = semiauto_change_image_callbacks(handles);

% initial draw
slider_callbacks_draw_image_slice(handles);

% a good order for doing error correction operations and the like
handles.layer_array = [handles.info.master_layer:my_sign(handles.info.top_layer-handles.info.bottom_layer):handles.info.top_layer ...
    handles.info.master_layer-my_sign(handles.info.top_layer-handles.info.bottom_layer):-my_sign(handles.info.top_layer-handles.info.bottom_layer):handles.info.bottom_layer];
handles.time_array = [handles.info.master_time:my_sign(handles.info.end_time-handles.info.start_time):handles.info.end_time ...
    handles.info.master_time-1:-my_sign(handles.info.end_time-handles.info.start_time):handles.info.start_time];

% make directories if they dont exist
[a b c] = mkdir(handles.tempsrc.parent);
[a b c] = mkdir(handles.tempsrc.bord);
% [a b c] = mkdir(handles.tempsrc.poly);

fprintf('Data Set "%s" loaded.\n', handles.data_set);
readyproc(handles, 'ready')