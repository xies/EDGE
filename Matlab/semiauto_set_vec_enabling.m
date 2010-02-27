function semiauto_set_vec_enabling(handles)

% if switching to automatic mode, unselect all
if get(handles.radiobutton_vec_manual, 'Value')
    vec_adjustments_visible(handles, 'on');

    % enable the vectorized cell adjustments
    names = fieldnames(handles);
    for i = 1:length(names)
        cont = strfind(names{i}, 'vec_');
        if ~isempty(cont) && cont(1) == 1
            set(handles.(names{i}), 'Enable', 'on')
        end
        cont = strfind(names{i}, 'button_refine_');
        if ~isempty(cont) && cont(1) == 1
            set(handles.(names{i}), 'Enable', 'on')
        end
    end

    set(handles.vec_activate_cell, 'Enable', 'off');
else
    handles.activeCell = [];
    handles.activeVertex = [];
    set(handles.vec_add_cell, 'Enable', 'off');
    set(handles.vec_remove_cell, 'Enable', 'off');
    set(handles.vec_move_vertex, 'Enable', 'off');
    set(handles.vec_remove_vertex, 'Enable', 'off');
    
    set(handles.vec_activate_cell, 'Enable', 'on');
    set(handles.vec_activate_cell, 'String', 'S+A+R stack');
end