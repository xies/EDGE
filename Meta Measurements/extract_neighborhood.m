function neigh = extract_neighborhood(data_set, cell_inds, max_order)

load(fullfile('/Users/mgelbart/Desktop/EDGE/DATA_GUI',data_set, 'embryo_data'));
% loads 'embryo4d'

if isempty(cell_inds)
    cell_inds = embryo4d.getCellGraph(embryo4d.masterTime, embryo4d.masterLayer).activeCellIndices;
end

neigh = cell(length(cell_inds), embryo4d.t, embryo4d.z, max_order);

for i = 1:length(cell_inds)
    j = 1;
    for time = embryo4d.startTime:sign(embryo4d.endTime-embryo4d.startTime):embryo4d.endTime
        k = 1;
        for layer = embryo4d.bottomLayer:sign(embryo4d.topLayer-embryo4d.bottomLayer):embryo4d.topLayer
            for l = 1:max_order
                neigh{i, j, k, l} = embryo4d.getCellGraph(time, layer).cellNeighborActiveIndices(i, l);
            end
            k = k + 1;
        end
        j = j + 1;
    end
end
    