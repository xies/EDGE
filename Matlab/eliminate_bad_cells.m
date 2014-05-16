% Clean up an image "input" and remove unwanted cells. First, perform
% erosions if necessary to remove junk on the outside of the image. Then, 
% remove cells with area smaller than min_cell_sz, which is given in
% pixel^2

function [dirty,clean] = eliminate_bad_cells(...
    input, min_cell_sz, num_erosions)

% remove junk on the outside
if num_erosions > 0
    for erode = 1:num_erosions
        input = fix_outer_cells_preserve_vertices(input);
    end
    % the below lines should not be here when the rest of the function
    % is uncommented, because it gets done at the bottom anyway!
%     dirty = bwmorph(input, 'skel', Inf);
%     clean = bwmorph(dirty, 'shrink', Inf);
%     clean = bwmorph(clean,  'clean',  Inf);
else
    dirty = input;
    clean = input;
end


% this part is removed for speed in favor of doing it in the CellGRaph
% constructor. however, it is unclear that the other method is always
% faster. if you have TONS of small cells, then it's better to get rid of
% them with this method that probably doesnt scaale that badly, as opposed
% to going through the CellGraph constructor with so many cells before
% finally removing the small ones at the last step.

regions = bwlabel(logical(1-input), 4);

% now it will not ignore any of the regions,
% since it sometimes confuses the borders as the background
regions = regions + 1;

tic
S = regionprops(regions, 'Area', 'Solidity');
toc
tic
cellAreas1 = [S.Area];
areas = [cellAreas1; 1:length(cellAreas1)]';


sol1 = [S.Solidity];
solidity = [sol1; 1:length(sol1)]';
toc
tic
[~, bordindex] = min(solidity(:,1));
toc

badcells = areas(areas(:,1) < min_cell_sz, 2);
badcells = badcells(:);

tic
% original = regions;
% make those into BORDERS and then skeletonize

[dirty,clean] = eliminate_region(regions,badcells(:),bordindex);
toc