function [data_allz data_topz data_btmz data_midz neighbor_indices] = ...
    extract_measurement(data_set, measurement, cell_inds, layers_from_top)

% load the data
data = [];
load(fullfile('/Users/mgelbart/Desktop/EDGE/DATA_GUI',data_set, 'Measurements', measurement));
% loads into "data", "name", "unit"

% load the cells of interest
cell_indices = [];
neighbor_indices = [];
if ~isempty(cell_inds)
    load(fullfile('/Users/mgelbart/Desktop/EDGE/DATA_OUTPUT',data_set,cell_inds));
else
    cell_indices = 1:size(data,3); 
end
% "loads into "cell_indices", "neighbor_indices"

% extract the data
data_allz = zeros(size(data, 1), size(data, 2), length(cell_indices));
data_topz = zeros(size(data, 1), length(cell_indices));
data_btmz = zeros(size(data, 1), length(cell_indices));
data_midz = zeros(size(data, 1), length(cell_indices));
for t = 1:size(data, 1)
    for i = 1:length(cell_indices)
        data_allz(t,:,i) = cell2mat(data(t,:,cell_indices(i)));
        
        % non-nan data
        good_data = find(~isnan(data_allz(t,:,i)));
        if length(good_data) <= layers_from_top
            data_topz(t, i) = NaN;
        else
            data_topz(t, i) = data_allz(t, good_data(end-layers_from_top),i); 
        end
        if length(good_data) < 1
            data_btmz(t, i) = NaN;
%             data_midz(t, i) = NaN;
        else
            data_btmz(t, i) = data_allz(t, good_data(1), i);
%             data_midz(t, i) = data_allz(t, good_data(ceil(length(good_data)/2)  ), i);
            % midz takes middle of non-NaN values, not middle of 
            % coordinate system.
        end
        data_midz(t, i) = data_allz(t, 9, i);  % temp1!
    end
end