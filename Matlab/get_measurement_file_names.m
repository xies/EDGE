function [measurechannels measures] = get_measurement_file_names(handles)
% finds the folder names of all the measurements that apply to the current
% data set

measure_path = fullfile(handles.program_dir, 'Measurements');

% channels = dir();
channels = handles.channelnames;
channels = ['Membranes'; channels(:)];

[measurechannels measures] = get_measurement_file_names_specific_ch(measure_path, channels);