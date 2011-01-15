% Finds the specific measurement files for a given data set. 
% Looks in the folder measure_path, and looks for the specified channels.
% Returns measures, and array of all measurement function names, and
% measurementchannels, an array of the same size that gives the
% corresponding channel for each of the measurement names. 

function [measurechannels measures] = get_measurement_file_names_specific_ch(measure_path, channels)

current_dir = pwd;

% first, we have the built in measurements
measures = cell(0);
measurechannels = cell(0);
% now we must check what is exported
% first, loop through each subdirector in the Measurements folder
% then loop through each file and get the file names
% 


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