function [clean cg] = semiauto_preprocess(handles, T, Z)

% if T == 51 && Z == 3
%     [clean cg] = semiauto_preprocess_active_contours(handles, T, Z);
%     return;
% end

% preprocess the given image in the semiautomatic image processing

cells = imread(handles.info.image_file(T, Z, handles.src.raw));
% cells = cells(:,:,1); % in case it's a 3 channel image

cells = double(cells);
% process into thin borders
lp = handles.info.bandpass_low / handles.info.microns_per_pixel;  % change to pixels
hp = handles.info.bandpass_high / handles.info.microns_per_pixel; % change to pixels
th = handles.info.preprocessing_threshold;

cellsi = get_membs_v3(cells, lp, hp, th);

% change the minimum cell size from microns^2 to pixels^2
minimum_cell_size = handles.info.minimum_cell_size / (handles.info.microns_per_pixel)^2;

% eliminate bad cells  
[dirty clean] = eliminate_bad_cells(...
    cellsi, minimum_cell_size, handles.info.number_of_erosions);

% get the cell properties (centroids, vertices, etc.)
[centroid_list regions] = find_centroids(clean);
[vertex_list] = find_vertices(dirty);

% the minimum distance between vertices before they get merged into one
% VERTEX_MERGE_DIST_THRESH = 3.0; %pixel
VERTEX_MERGE_DIST_THRESH = handles.info.refine_min_edge_length / handles.info.microns_per_pixel;
% convert from microns to pixels
% note that the min_edge_length is used here and when doing automatic split
% edge -- in both cases, one cannot make edges smaller than this
% also, note that the VERTEX_MERGE_DIST_THRESH may not end up being exactly
% the same as a minimum edge length, but it's close enough

% the minimum angle formed by vertices
VERTEX_MIN_ANGLE_THRESH = degtorad(handles.info.refine_min_angle);
% convert it to radians
% note this is also used for refining edges, but we enforce it here
% by removing vertices that cause a violation

% create the CellGraph object
cg = CellGraph(regions, centroid_list, vertex_list, T, Z, ...
    VERTEX_MERGE_DIST_THRESH, VERTEX_MIN_ANGLE_THRESH);
