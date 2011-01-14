function handles = slider_callbacks_draw_image_slice_dots_semiauto(handles)
% the draw_image_slice function is now split into two. this one deals with
% when you click on cells or vertices-- it just plots the dots as
% necessary, without replotting the rest of handles.axes1.


axes(handles.axes1);
hold on


try  % delete previous centroids
    delete(handles.cents_handle);
end
try  % delete previous vertices
    delete(handles.verts_handle);
end

% draw centroids of active Cells
if ~isempty(handles.activeCell)
    % make a javaarray
%     activecells = javaArray('Cell', length(handles.activeCell));
%     for i = 1:length(handles.activeCell)
%         activecells(i) = handles.activeCell{i};
%     end

    [T Z] = getTZ(handles);
    activecells = handles.embryo.getCells(handles.activeCell, T, Z);
    cents = Cell.centroidStack(activecells);
    handles.cents_handle = plot(cents(:,2), cents(:,1), '.r');  
    
    verts = Cell.vertexCoords(activecells);
    handles.verts_handle = plot(verts(:,2), verts(:,1), '.b');  %'MarkerSize', 15
    
    % plot the vertices of these cells
%     for i = 1:length(activecells)
%         verts = Vertex.coords(activecells(i).vertices);
%         handles.verts_handle = plot(verts(:,2), verts(:, 1), '.r');
%     end
end



% for semiauto, draw vertices
if isfield(handles, 'activeVertex')
    if ~isempty(handles.activeVertex)
        % make a javaarray
        activeverts = javaArray('Vertex', length(handles.activeVertex));
        for i = 1:length(handles.activeVertex)
            activeverts(i) = handles.activeVertex{i};
        end
        verts = Vertex.coords(activeverts);
        handles.verts_handle = plot(verts(:,2), verts(:, 1), '.b');
    end
end
