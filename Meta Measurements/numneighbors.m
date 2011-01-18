clear all; close all; clc;

%%
% data_set = 'Histone032';
data_set = '090309 MembCherry HistoneGFP';
cell_inds = []; % all cells
measurement = 'Membranes--vertices--# of neighbors';
layers_from_top = 3;

%%
% '1-15-2011, 3;29 PM'

[data_all data_apical data_basal data_middle] = ...
    extract_measurement(data_set, measurement, cell_inds, layers_from_top);

%% take distribution over cells
out = [];
%% 
out(:,1)=neighbor_results(data_apical, 'apical');
%%
out(:,2)=neighbor_results(data_basal,  'basal');
%%
out(:,3)=neighbor_results(data_middle, 'middle');


%% matthias wants the image thing. OK
plot_neighbor_image(diff(data_apical), 5);
%%
plot_neighbor_image(data_middle, 16);
%%
plot_neighbor_image(diff(data_basal),  7);