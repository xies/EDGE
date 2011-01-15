function plot_measurement_vs_AP_axis_at_apical_end(data_set, measurement, cell_inds)


[data_allz data_topz] = extract_measurement(data_set, measurement, cell_inds);

% also need to load in the coordinate to sort by AP. this is why we wanted
% them to already be oriented...! but that's ok, sort by x-coord is still
% basically correct.
centx_name = 'Membranes--basic_2d--Centroid-x';

[centx_allz centx_topz] = extract_measurement(data_set, centx_name, cell_inds);



to_sort = [centx_topz(:) data_topz(:)];
sorted = sortrows(to_sort);

figure;
plot(sorted(:,1),sorted(:,2));
xlabel('A-P')
ylabel(['Apical ' measurement]);
title(['Apical ' measurement ' for cells near midline along A-P axis']);