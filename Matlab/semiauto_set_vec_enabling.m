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

    set(handles.vec_activate_cell, 'String', 'Activate cell');
    set(handles.vec_activate_cell, 'Enable', 'on');

%     set(handles.vec_activate_cell, 'Enable', 'off');    
%     [T Z] = getTZ(handles);
%     for i = 1:length(handles.activeCell)
%         if handles.embryo.isTrackingCandidate(handles.activeCell(i), T, Z)
%             set(handles.vec_activate_cell, 'Enable', 'on');
%         end
%     end
    
else
    handles.activeCell = [];
    handles.activeVertex = [];
    set(handles.vec_add_cell, 'Enable', 'off');
    set(handles.vec_remove_cell, 'Enable', 'off');
    set(handles.vec_move_vertex, 'Enable', 'off');
    set(handles.vec_remove_vertex, 'Enable', 'off');
    set(handles.vec_puncture_cell, 'Enable', 'off');
    
    set(handles.vec_activate_cell, 'Enable', 'on');
    set(handles.vec_activate_cell, 'String', 'R+S+A edge');
end

% this is screwy-- should only add cells by adding edges
% it doesn't work because if i cut into a cell to add a new one, the
% vertices of that existing cell don't get updated. basically, this
% function does not take into account the effects it has on existing cells,
% but just throws a new one in the mix. AddEdge does this properly with no
% problems
% set(handles.vec_add_cell, 'Enable', 'off');