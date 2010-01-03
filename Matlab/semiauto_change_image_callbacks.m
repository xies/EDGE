function handles = semiauto_change_image_callbacks(handles)
% resets the temp CG file every time you change the image, and other things
% should be called whenever you change what image you're looking at in
% semiauto

[T Z] = getTZ(handles);
% handles.tempcg = handles.embryo.getCellGraph(T, Z);

if isempty(handles.embryo.getCellGraph(T, Z))
    
    % disable the vectorized cell adjustments
    names = fieldnames(handles);
    for i = 1:length(names)
        cont = strfind(names{i}, 'vec_'); % if it's one of these buttons
        if ~isempty(cont) && cont(1) == 1
            set(handles.(names{i}), 'Enable', 'off')
        end
        cont = strfind(names{i}, 'button_refine_');
        if ~isempty(cont) && cont(1) == 1
            set(handles.(names{i}), 'Enable', 'off')
        end
    end
else    
    
    semiauto_set_vec_enabling(handles);
%     % enable the vectorized cell adjustments
%     names = fieldnames(handles);
%     for i = 1:length(names)
%         cont = strfind(names{i}, 'vec_');
%         if ~isempty(cont) && cont(1) == 1
%             set(handles.(names{i}), 'Enable', 'on')
%         end
%         cont = strfind(names{i}, 'button_refine_');
%         if ~isempty(cont) && cont(1) == 1
%             set(handles.(names{i}), 'Enable', 'on')
%         end
%     end
end   

handles.activeCell = [];
handles.activeVertex = [];


% enable/disable as needed
if handles.embryo.isTracked
    set(handles.button_export, 'Enable', 'on');
%     set(handles.button_refine_applyall, 'Enable', 'on');
    set(handles.cbox_inactive, 'Enable', 'on');
else
    set(handles.button_export, 'Enable', 'off');
%     set(handles.button_refine_applyall, 'Enable', 'off');
    set(handles.cbox_inactive, 'Enable', 'off');
end


