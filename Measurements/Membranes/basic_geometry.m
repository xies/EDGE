function [data names units] = basic_geometry(embryo, getMembranes, t, z, c, dx, dz, dt, other)
% computes the basic geometrical properties at (t, z) for Cell c given the
% Embryo4D and the relevant resolutions

names{1} = 'Cell length';
names{2} = 'Cell volume';

units{1} = 'microns';
units{2} = 'microns^3';

% % just do this for one layer to save time, since these properties are a
% % function of the whole cell stack, not an individual layer
% if z ~= embryo.masterLayer
%     data{1} = NaN;
%     data{2} = NaN;
%     return;
% end

if ~embryo.anyTracked(c, t)
    data{1} = NaN;
    data{2} = NaN;
    return;
end

%% compute the cell length
centroids = Cell.centroidStack(embryo.getCellStack(c, t)) * dx;
y_values = centroids(:, 1);
x_values = centroids(:, 2);
z_values = (1:length(x_values)) * dz;
z_values = z_values(:);

% remove NaN values
for i = length(x_values):-1:1
    if isnan(x_values(i))
        x_values(i) = [];
        y_values(i) = [];
        z_values(i) = [];
    end
end
cell_length = sum(sqrt(diff(x_values).^2 + ...
                           diff(y_values).^2 + ...
                           diff(z_values).^2));

%% compute the cell volume 
% here we just compute each volume element as area * dz. we might
% be able to do a smarter estimation using the orthogonal area, but this is
% ok for now.

% add up the volume elements from each layer. if the current
% layer is not tracked (NaN), then use the last volume element that you
% found. move from the bottom to the top. do not include an element from
% the top, because there is nothing "above" it
cell_volume = 0;
dir = sign(embryo.highestTracked(c, t) - embryo.lowestTracked(c, t));
for i = embryo.lowestTracked(c, t) : dir : embryo.highestTracked(c, t) - dir
   if ~isempty(embryo.getCell(c, t, i));
       last_volume_element = (embryo.getCell(c, t, i).area * dx^2) * dz;
   end
   cell_volume = cell_volume + last_volume_element;
end


%% places the data in the cell array

data{1} = cell_length;
data{2} = cell_volume;

