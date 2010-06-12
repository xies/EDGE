function [cents regions] = find_centroids(bords)

% cannot use bwconncomp here because we need to pass the regions graph to
% the CellGraph constructor.
regions = bwlabel(logical(1-bords), 4);
cent = regionprops(regions, 'Centroid');
cents = [cent.Centroid];
cents = reshape(cents, [2 length(cents)/2]).';
cents = round(cents);

% flip x and y
cents = fliplr(cents);

%gets rid of the background
% cents = removerows(cents, 1); 
cents(1,:)=[];

%fix region array
regions = regions - 1;
