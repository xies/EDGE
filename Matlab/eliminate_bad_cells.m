function [dirty clean] = eliminate_bad_cells(...
    cellsA, min_cell_sz, NUM_OF_EROSIONS)

% note: **min_cell_sz is given in pixel^2

MINIMUM_SIZE = min_cell_sz; %pixel^2

% this is under construction. i started using
% fix_outer_cells_preserve_vertices and i think it is a very good idea
% because cells get punctured elsewhere anyway so no need to destroy the
% vertices. 
% another idea would be to put the fix_oter_cells parts BEFORe the
% detruction of overly small cells, and/or to use a cell size threshhold
% based on the average cell size instead of a fixel pixel threshhold that
% only depends on time
% figure, imshow(cellsA);
% keyboard;

for erode = 1:NUM_OF_EROSIONS
% remove junk on the outside %****
    cellsA = fix_outer_cells_preserve_vertices(cellsA);
end
cellsC = cellsA;

% figure, imshow(cellsC)
% keyboard
%x = fix_outer_cells(x);
%bords = fix_outer_cells(x);  % kills the outermost layer

% trying to use a lookup table instead of fix outer cells
%{
regions = bwlabel(logical(1-x), 4);

S = regionprops(regions,'Area');
cellAreas1 = [S.Area];
areas = [cellAreas1; 1:length(cellAreas1)]';
% make bgindex equal zero
[dummy bgindex] = max(areas(:,1));

regions(regions == bgindex) = (1/2);

% instead of::: bords = fix_out_cells(x);
f = @(x) (sum( x(:) == (1/2) ) == 0);
lut = makelut(f, 3);
regions = applylut(x, lut);

bords = ~regions;
%}




regions = bwlabel(logical(1-cellsC), 4);
% makes the borders be the bacground

% now it will not ignore any of the regions,
% since it sometimes confuses the borders as the background
regions = regions + 1;

% max(regions(:))
% 
% figure, imshow(regions == 1);
% title 'equals1'
% figure, imshow(regions == 0);
% title 'equals0'


S = regionprops(regions, 'Area', 'Solidity');
cellAreas1 = [S.Area];
areas = [cellAreas1; 1:length(cellAreas1)]';

% remove the background from the areas list...
[dummy bgindex] = max(areas(:,1));
%areas = removerows(areas, bgindex);




sol1 = [S.Solidity];
solidity = [sol1; 1:length(sol1)]';
[dummy bordindex] = min(solidity(:,1));




badcells2 = areas(areas(:,1) < MINIMUM_SIZE, 2);
badcells2 = badcells2(:);

%badcells = [badcells1; badcells2];
badcells = badcells2;


% make those into BORDERS and then skeletonize
for i = 1:length(badcells)
    regions(regions == badcells(i)) = bordindex;
end

cellsD = regions == bordindex;


%{
figure, imshow(finalbords)
title 'after fix outer cells'
% now, round 2, get rid of cells that are a standard deviation area less
% than others
regions = bwlabel(logical(1-finalbords), 4);
regions = regions + 1;

S = regionprops(regions,'Area','Solidity');
cellAreas1 = [S.Area];
areas = [cellAreas1; 1:length(cellAreas1)]';

% remove the background and borders from the areas list...
[dummy bgindex] = max(areas(:,1));

sol1 = [S.Solidity];
solidity = [sol1; 1:length(sol1)]';
[dummy bordindex] = min(solidity(:,1));

areas = removerows(areas, [bgindex bordindex]);
AREA_THRESH_DEV = 3.0;
AREA_THRESH = median(areas(:,1)) - AREA_THRESH_DEV * std(areas(:,1));
sort(areas(:,1))
median(areas(:,1))
AREA_THRESH

badcells1 = areas(areas(:,1) < AREA_THRESH, 2);
badcells1 = badcells1(:);
for i = 1:length(badcells1)
    regions(regions == badcells1(i)) = bgindex;  % set to BG instead of BORD becauase these should touch outside
end
newbords = regions == bordindex;
finalbords = newbords;

figure, imshow(finalbords)
%}


cellsE = bwmorph(cellsD, 'skel', Inf);

clean = bwmorph(cellsE, 'shrink', Inf);
clean = bwmorph(clean,  'clean',  Inf);

dirty = cellsE;

% need to get rid of some islands on "dirty"
% do this by finding the region of biggest area (aka the main borders)
% and only keeping that. this should solve some problems. note that this
% is much easier than using sort, which I think I was doing unnecessesarily
% before...

% NEW: this is a problem when are are multiple regions
% of cells, such as when you go very deep in the fixed
% data sets!!!!!! therefore we should not do it this way
% but rather by removing islands??
% % % % dirty = bwlabel(dirty);
% % % % props = regionprops(dirty, 'Area');
% % % % props = [props.Area].';
% % % % [val ind] = max(props);
% % % % dirty = dirty == ind;

% figure, imshow(finalbords)
% title 'after skel'


%{
% ADDITIONAL CODE:
% we want to get rid of things sticking out in the middle of cells, while
% still keeping the things sticking out at the edges because they give the
% vertices of edge cells

mask = bwlabel(logical(1-finalbords), 4) > 1;
mask = bwmorph(mask, 'dilate');

inner = bwmorph(finalbords.*mask,'shrink', Inf);
inner = bwmorph(inner,'clean');

%patch them together so that the shrink and clean only act on the inside
finalbords = finalbords .* (1-mask) + inner.*mask;
%}
