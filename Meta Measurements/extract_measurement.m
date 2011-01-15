function [data_allz data_topz] = extract_measurement(data_set, measurement, cell_inds, layers_from_top)

% load the data
load(fullfile('/Users/mgelbart/Desktop/EDGE/DATA_GUI',data_set, 'Measurements', measurement));
% loads into "data", "name", "unit"

% load the cells of interest
if ~isempty(cell_inds)
    load(fullfile('/Users/mgelbart/Desktop/EDGE/DATA_OUTPUT',data_set,cell_inds));
else
    cell_indices = 1:size(data,3); 
end
% "loads into "cell_indices"

% extract the data
data_allz = zeros(size(data, 1), size(data, 2), length(cell_indices));
data_topz = zeros(size(data, 1), length(cell_indices));
for t = 1:size(data, 1)
    for i = 1:length(cell_indices)
        data_allz(t,:,i) = cell2mat(data(t,:,cell_indices(i)));
        apical_data = find(~isnan(data_allz(t,:,i)));
        if length(apical_data) <= layers_from_top
            data_topz(t, i) = NaN;
        else
            data_topz(t, i) = data_allz(t, apical_data(end-layers_from_top),i); 
        end
    end
end