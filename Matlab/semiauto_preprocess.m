function [clean cg] = semiauto_preprocess(handles, T, Z)

% preprocess the given image in the semiautomatic image processing

cells = imread(handles.info.image_file(T, Z, handles.src.raw));
cells = cells(:,:,1); % in case it's a 3 channel image

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
VERTEX_MERGE_DIST_THRESH = 3.0; %pixel

% create the CellGraph object
cg = CellGraph(regions, centroid_list, vertex_list, T, Z, VERTEX_MERGE_DIST_THRESH);
