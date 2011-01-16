clear all; close all; clc;


data_set = 'Histone032';
cell_inds = []; % all cells
measurement = 'Membranes--vertices--# of neighbors';
layers_from_top = 3;


% '1-15-2011, 3;29 PM'

[data_all data_apical] = extract_measurement(data_set, measurement, cell_inds, layers_from_top);

%%


%average over cells, ignoring NaNs
meannn = my_mean(data_apical);

% plot # neighbors vs time
figure;
plot(meannn);

figure
plot(data_apical)



