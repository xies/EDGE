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

% handles.activeCell = [];
% some of the indices that are slected might not be in the new image, so
% unselect those just for simplicity. also, i can only do this for _active_
% cells because inactive cells might have the same index between images but
% have nothing to do with each other
if handles.embryo.isTracked
    handles.activeCell = intersect(handles.activeCell, handles.embryo.getCellGraph(T, Z).activeCellIndices);
    if length(handles.activeCell) == 1  % fix the text
            set(handles.cell_text, 'String', num2str(handles.activeCell(1))); 
        else
            set(handles.cell_text, 'String', '-'); 
    end
else
    handles.activeCell = [];
    set(handles.cell_text, 'String', '-'); 
end
handles.activeVertex = [];


% enable/disable buttons and checkboxes as needed
if handles.embryo.isTracked
    set(handles.button_export, 'Enable', 'on');
    set(handles.cbox_inactive, 'Enable', 'on');
else
    set(handles.button_export, 'Enable', 'off');
    set(handles.cbox_inactive, 'Enable', 'off');
end
if exist(handles.info.image_file(T, Z, handles.tempsrc.bord), 'file')
    % if they have already processed that image, the file will exist
    set(handles.cbox_bord, 'Enable', 'on');
else
    set(handles.cbox_bord, 'Enable', 'off');
end
if ~isempty(handles.embryo.getCellGraph(T, Z))
    % if there is a polygon, should be same as above but i check separately 
    set(handles.cbox_poly, 'Enable', 'on');
else
    set(handles.cbox_poly, 'Enable', 'off');
end


