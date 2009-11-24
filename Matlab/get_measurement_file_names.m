function [measurechannels measures] = get_measurement_file_names(handles)
% finds the folenames of all the measurements in the measurements directory
% and all subfolders. removes the .m from the end of each one


current_dir = pwd;

% first, we have the built in measurements
measures = cell(0);
measurechannels = cell(0);
% now we must check what is exported
% first, loop through each subdirector in the Measurements folder
% then loop through each file and get the file names
% 
measure_path = fullfile(handles.program_dir, 'Measurements');

% channels = dir(measure_path);
channels = handles.channelnames;
channels = ['Membranes'; channels(:)];

for i = 1:length(channels)
    cd(fullfile(measure_path, channels{i}));  % go in that directory
    functions = dir();
    for j = 1:length(functions)
        if functions(j).isdir
            continue;
        end
        fname = functions(j).name;
        if length(fname) <= 2
            continue;
        end
        if ~strcmp(fname(end-1:end), '.m')  % make sure it ends with '.m'
            continue;
        end
            
        function_name = functions(j).name;
        function_name = function_name(1:end-2);  % remove the .m from the end

        measures = [measures function_name];
        measurechannels = [measurechannels channels{i}];
    end
end

cd(current_dir);