function handles = handles_set_program_dir(handles)

program_dir = mfilename('fullpath');
% remove the 'Matlab' subdirectory from the end of the path
matlab_in_string = strfind(program_dir, 'Matlab');
matlab_in_string = matlab_in_string(end); % just in case 'Matlab' appears upstream in the path
program_dir = program_dir(1:matlab_in_string-2);
handles.program_dir = program_dir;  % the directory of the whole program

% could do this better by doing cd('..') then program_dir = dir, then
% cd('Matlab')
