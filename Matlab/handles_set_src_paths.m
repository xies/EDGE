function handles = handles_set_src_paths(handles)

source = fullfile(handles.program_dir, 'DATA_GUI', handles.data_set);

handles.src.parent = source;        % the subdirectory for this data set
handles.src.membranes = fullfile(source, 'Membranes');
handles.src.raw = fullfile(source, 'Membranes', 'Raw');  % subsub dirctory: raw
handles.src.bord = fullfile(source, 'Membranes', 'Processed');
handles.src.measurements = fullfile(source, 'Measurements');