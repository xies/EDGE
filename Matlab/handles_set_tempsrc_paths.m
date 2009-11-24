function handles = handles_set_tempsrc_paths(handles)

% this below is for the semiautomatic processing temporary files
source = fullfile(handles.program_dir, 'DATA_SEMIAUTO', handles.data_set);

handles.tempsrc.parent = source;
handles.tempsrc.bord = fullfile(source, 'Processed membranes');
% handles.tempsrc.poly = fullfile(source, 'cg Objects');
