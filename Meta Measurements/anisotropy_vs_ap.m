clear all; close all; clc;

% plot some stuff for the paper.
% using fixed data set c1_20 for now, plot anisotropy

%%
% data_set = 'c1_20';  
data_set = 'c2_6';
% cell_inds = '1-14-2011, 1;12 PM';
cell_inds = '1-14-2011, 2;55 PM';
% cell_inds = []; % all cells
%%

measurement = 'Membranes--ellipse_properties--Anisotropy';
plot_measurement_vs_AP_axis_at_apical_end(data_set, measurement, cell_inds)
plot_measurement_vs_DV_axis(data_set, measurement, cell_inds)

%%
measurement = 'Membranes--vertices--# of vertices';
plot_measurement_vs_AP_axis_at_apical_end(data_set, measurement, cell_inds)
plot_measurement_vs_DV_axis(data_set, measurement, cell_inds)

%%
measurement = 'Membranes--vertices--# of neighbors';
plot_measurement_vs_AP_axis_at_apical_end(data_set, measurement, cell_inds)
plot_measurement_vs_DV_axis(data_set, measurement, cell_inds)

%%

measurement = 'Membranes--basic_2d--Area';
[data_allz data_topz] = extract_measurement(data_set, measurement, cell_inds);
centx_name = 'Membranes--basic_2d--Centroid-x';
[centx_allz centx_topz] = extract_measurement(data_set, centx_name, cell_inds);

% calculate amount of apical constriction
apicalc = zeros(size(data_allz, 2),1);
for i = 1:size(data_allz, 2)
    x=data_allz(:,i);
    x=x(~isnan(x));
    apicalc(i) = 1/(sum(x(1:5))/sum(x(end-4:end)));
end

to_sort = [centx_topz(:) apicalc(:)];
sorted = sortrows(to_sort);

figure;
plot(sorted(:,1),sorted(:,2));
xlabel('A-P')
ylabel(['apical constriction']);
