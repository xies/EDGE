% Finds the filenames of all the measurements that apply to the current
% data set by searching in the EDGE/Measurements folder.

function [measurechannels measures] = get_measurement_file_names(handles)

measure_path = fullfile(handles.program_dir, 'Measurements');

channels = handles.channelnames;
channels = ['Membranes'; channels(:)];

[measurechannels measures] = get_measurement_file_names_specific_ch(measure_path, channels);