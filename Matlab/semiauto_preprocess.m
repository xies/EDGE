function [clean cg] = semiauto_preprocess(handles, T, Z)

% if T == 51 && Z == 3
%     [clean cg] = semiauto_preprocess_active_contours(handles, T, Z);
%     return;
% end

t = zeros(1,6);

% preprocess the given image in the semiautomatic image processing
tic
cells = imread(handles.info.image_file(T, Z, handles.src.raw));
t(1)=toc;
% cells = cells(:,:,1); % in case it's a 3 channel image


% a special clause for totally blank images (these are "fillers")
if ~any(cells(:))
   CellGraph(cells, [], [], T, Z, NaN, NaN);
end



cells = double(cells);
% process into thin borders
lp = handles.info.bandpass_low / handles.info.microns_per_pixel;  % change to pixels
hp = handles.info.bandpass_high / handles.info.microns_per_pixel; % change to pixels
th = handles.info.preprocessing_threshold;

tic
cellsi = get_membs_disperse(handles.info.image_file(T, Z, handles.src.raw), th, handles.info.Xs, handles.info.Ys);
t(2)=toc;

% change the minimum cell size from microns^2 to pixels^2
minimum_cell_size = handles.info.minimum_cell_size / (handles.info.microns_per_pixel)^2;

% tic
% eliminate bad cells  
[dirty clean] = eliminate_bad_cells(...
    cellsi, minimum_cell_size, handles.info.number_of_erosions);
% t(3)=toc;

% get the cell properties (centroids, vertices, etc.)
tic
[centroid_list regions] = find_centroids(clean);
t(4)=toc;
tic
[vertex_list] = find_vertices(dirty);
t(5)=toc;

% the minimum distance between vertices before they get merged into one
% VERTEX_MERGE_DIST_THRESH = 3.0; %pixel
VERTEX_MERGE_DIST_THRESH = handles.info.refine_min_edge_length / handles.info.microns_per_pixel;
% convert from microns to pixels
% note that the min_edge_length is used here and when doing automatic split
% edge -- in both cases, one cannot make edges smaller than this
% also, note that the VERTEX_MERGE_DIST_THRESH may not end up being exactly
% the same as a minimum edge length, but it's close enough
 
%yes, using it for both is bad - it should be higher for the error
%correction than here. same thing as below with the angle thresh
% i know hardcoding is bad, but ... oh well ;)
VERTEX_MERGE_DIST_THRESH = 2; % pixels?!?
% ok, now i give a negative one because we use 'branchpoints' to find
% vertices in find_vertices.m
% this should be adjustable!!



% the minimum angle formed by vertices
VERTEX_MIN_ANGLE_THRESH = degtorad(handles.info.refine_min_angle);
% convert it to radians
% note this is also used for refining edges, but we enforce it here
% by removing vertices that cause a violation
% actually, this shouldn't be the same as in refining edges...!
% so we just use something else, a default value of 15 degrees
VERTEX_MIN_ANGLE_THRESH = degtorad(15);

% create the CellGraph object
tic
cg = CellGraph(regions, centroid_list, vertex_list, T, Z, ...
    VERTEX_MERGE_DIST_THRESH, VERTEX_MIN_ANGLE_THRESH, minimum_cell_size);
t(6)=toc;
