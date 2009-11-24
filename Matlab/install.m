% Install script for EDGE
% Michael Gelbart // 2009
clear all; clc;

% ***user-defined *** %
JAVA_NAME = fullfile('java1', 'bin');
MATLAB_NAME = 'Matlab';
MEASUREMENTS_NAME = 'Measurements';
JAVA_MEMORY_STRING = '-Xmx512m';
% ******************* %


% get the path to this file
program_dir = mfilename('fullpath');
% remove the 'Matlab' subdirectory from the end of the path
matlab_in_string = strfind(program_dir, MATLAB_NAME);
matlab_in_string = matlab_in_string(end); % just in case 'Matlab' appears upstream in the path
program_dir = program_dir(1:matlab_in_string-1);
% path to the java
java_path = strcat(program_dir, JAVA_NAME);


%% edit the java.opts file
root_dir = matlabroot;  % the Matlab root directory
% opts_directory_location = strcat(root_dir, '/bin/', computer('arch'), '/java.opts');
opts_directory_location = fullfile(root_dir, 'bin', computer('arch'), 'java.opts');
fid = fopen(opts_directory_location, 'a');
if fid < 0
    disp(strcat('install.m: Cannot find file ', opts_directory_location));
    disp('Installation aborted.');
    return;    
end
fprintf(fid, '\n%s', JAVA_MEMORY_STRING);
fclose(fid);

%% edit the classpath file
% classpath_name = strcat(fullfile(matlabroot,'toolbox', 'local'), '/classpath.txt');
classpath_name = fullfile(matlabroot, 'toolbox', 'local', 'classpath.txt');
fid = fopen(classpath_name, 'a');
if fid < 0
    disp(strcat('install.m: Cannot find file ', classpath_name));
    disp('Installation aborted.');
    return;
end
fprintf(fid, '\n%s', java_path);
fclose(fid);

%% add the Measurements and Matlab folders to the path
% matlab_path = strcat(program_dir, MATLAB_NAME, '/');
matlab_path = fullfile(program_dir, MATLAB_NAME);
% measurements_path = strcat(program_dir, MEASUREMENTS_NAME, '/');
measurements_path = fullfile(program_dir, MEASUREMENTS_NAME);
addpath(matlab_path, measurements_path);
savepath;


%% message to user
msgbox('Installation completed. You must restart Matlab before continuing.', ...
    'Installation successful');

