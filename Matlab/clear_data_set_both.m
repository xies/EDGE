function handles = clear_data_set_both(handles)

% find all the other channels that are not Membranes
handles.channelnames = get_folder_names(handles.src.parent);
handles.channelnames(strcmp(handles.channelnames, 'Membranes')) = [];



% if only 1 other channel, set the name of that checkbox to be that
if length(handles.channelnames) == 1
    set(handles.cbox_other, 'String', handles.channelnames{1});
    set(handles.cbox_other, 'Enable', 'on'); 
elseif length(handles.channelnames) > 1
    set(handles.cbox_other, 'String', 'Other channels');
    set(handles.cbox_other, 'Enable', 'on'); 
else
    set(handles.cbox_other, 'Enable', 'off'); 
end


handles = handles_set_channelsrc_paths(handles);

% all the channels that have ever been exported
% used in slider_callbacks_draw_image_slice
handles.all_channelnames = get_folder_names(fullfile('..', 'Measurements'));

handles.activeChannels = [];

% is the data set fixed?
handles.fixed = isnan(handles.info.seconds_per_frame);


% the name of the function that gives the filename for an image
cd(handles.src.membranes);
handles.info.image_file = @image_filename;
% handles.info.image_file.raw = @(t, z) double(imread(image_filename(t, z, handles.src.raw)));
% handles.info.image_file.bord= @(t, z) double(imread(image_filename(t, z, handles.src.bord)));
handles.info.channel_image_file = cell(length(handles.channelnames));
for i = 1:length(handles.channelnames)
    cd(handles.src.channelsrc{i});
    handles.info.channel_image_file{i} = @image_filename;
%     handles.info.channel_image_file{i} = @(t, z) double(imread(image_filename(t, z, handles.src.channelsrc{channelnum})));
end
cd(fullfile(handles.program_dir, 'Matlab'));

% clear current axes
cla(handles.axes1);

% set the size of the axes
set(handles.axes1, 'XLim', [1 handles.info.Xs]);
set(handles.axes1, 'YLim', [1 handles.info.Ys]);

% set the "data info" panel
names = fieldnames(handles.info);
for i = 1:length(names)
    name = names{i};
    button_name = strcat('info_text_', name);
    if isfield(handles, button_name)
        set(handles.(button_name), 'String', my_num2str(handles.info.(name)));
    end
end

% initialize the checkboxes
set(handles.cbox_raw,  'Value', 1);
set(handles.cbox_bord, 'Value', 0);
set(handles.cbox_poly, 'Value', 1);
set(handles.cbox_other,'Value', 0);
%set(handles.cbox_inactive  -- only for semiauto?

