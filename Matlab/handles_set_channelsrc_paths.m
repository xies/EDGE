function handles = handles_set_channelsrc_paths(handles)

channelsrc = cell(length(handles.channelnames));
for i = 1:length(handles.channelnames)
    channelsrc{i} = fullfile(handles.src.parent, handles.channelnames{i});
end
handles.src.channelsrc = channelsrc;