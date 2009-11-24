function handles = vec_adjustments_visible(handles, onoff)

  
names = fieldnames(handles);
for i = 1:length(names)
    cont = strfind(names{i}, 'vec_');
    if ~isempty(cont) && cont(1) == 1
        set(handles.(names{i}), 'Visible', onoff)
    end
end

if strcmp(onoff, 'on')
    set(handles.button_vec_cancel, 'Visible', 'off');
else
    set(handles.button_vec_cancel, 'Visible', 'on');
end

if strcmp(handles.activeAdjustment, 'add_vertex')
    set(handles.button_vec_cancel, 'String', 'Done');
else
    set(handles.button_vec_cancel, 'String', 'Cancel');
end