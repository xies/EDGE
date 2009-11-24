function handles = initialize_tz_sliders(handles)

% set the text for z
start_layer = abs(handles.info.master_layer - handles.info.bottom_layer);
set(handles.z_slider, 'Value', start_layer);
set(handles.z_text, 'String', num2str(handles.info.master_layer));

% set the step values
z_max = abs(handles.info.bottom_layer - handles.info.top_layer);
set(handles.z_slider, 'Max', z_max);
set(handles.z_slider, 'Min', 0);
set(handles.z_slider, 'SliderStep', ... 
    [1/z_max; 1/z_max]);

% turn temporal slider off if fixed data set
if ~handles.fixed
    % set the text for t
    start_time = abs(handles.info.master_time - handles.info.start_time);
    set(handles.t_slider, 'Value', start_time);
    set(handles.t_text, 'String', num2str(handles.info.master_time));
    
    % times series data set
    set(handles.t_slider, 'Visible', 'on');
    set(handles.t_text_label, 'Visible', 'on');
    set(handles.t_text, 'Visible', 'on');
    
    t_max = handles.info.end_time - handles.info.start_time;
    set(handles.t_slider, 'Max', t_max);
    set(handles.t_slider, 'Min', 0);
    set(handles.t_slider, 'SliderStep', ...
        [1/t_max; 1/t_max]);
else
    % fixed data set
    set(handles.t_slider, 'Visible', 'off');
    set(handles.t_text_label, 'Visible', 'off');
    set(handles.t_text, 'Visible', 'off');
end

% need this because they are by default disabled when running a 
% new data set for the first time
set(handles.z_slider, 'Enable', 'on');
set(handles.t_slider, 'Enable', 'on');