function handles = go_to_image(handles, T, Z)

% called by update_embryo~~
set(handles.t_text, 'String', num2str(T));
set(handles.z_text, 'String', num2str(Z));

t_slider_value = abs(T - handles.info.start_time);
z_slider_value = abs(Z - handles.info.bottom_layer);
set(handles.t_slider, 'Value', t_slider_value);
set(handles.z_slider, 'Value', z_slider_value);


handles = semiauto_change_image_callbacks(handles);

% redraw
handles = slider_callbacks_draw_image_slice(handles);
