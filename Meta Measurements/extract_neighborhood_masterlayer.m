function neigh = extract_neighborhood_masterlayer(data_set, cell_inds, max_order)

load(fullfile('/Users/mgelbart/Desktop/EDGE/DATA_GUI',data_set, 'embryo_data'));
% loads 'embryo4d'

if isempty(cell_inds)
    cell_inds = embryo4d.getCellGraph(embryo4d.masterTime,embryo4d.masterLayer).activeCellIndices;
end

neigh = cell(length(cell_inds), embryo4d.t, max_order);

for i = 1:length(cell_inds)
    i
    j = 1;
    for time = embryo4d.startTime:sign(embryo4d.endTime-embryo4d.startTime):embryo4d.endTime
        for l = 1:max_order
            neigh{i, j, l} = embryo4d.getCellGraph(...
                time, embryo4d.masterLayer).cellNeighborActiveIndices(i, l);
        end
        j = j + 1;
    end
end
    