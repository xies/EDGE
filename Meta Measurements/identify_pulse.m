clear all; close all; clc;

% measurement that identifies apical constriction pulses

data_set = 'Histone032';
cell_inds = []; % all cells
measurement = 'Membranes--basic_2d--Area';
layers_from_top = 3;


% '1-15-2011, 3;29 PM'

[data_all data_apical] = extract_measurement(data_set, measurement, cell_inds, layers_from_top);

%%

% looking at the apical data only
% first, compute the rate of change of area
rate = diff(data_apical, 1);  % over time

%%

rate_nonan = rate(~isnan(rate));

% then, identify pulses by thresholding this rate at 1 st above meant
pulses = rate > mean(rate_nonan) + 2*std(rate_nonan);

%%
%average over cells, ignoring NaNs
meanrate = my_mean(rate);

% plot average constriction rate vs time
figure;
plot(meanrate);


% plot pulses vs time (just a thresholded version)
figure; 
imagesc(pulses);
xlabel('cell index')
ylabel('time');

% want to do a correlation based on neighbors
% where is the neighbor info stored? i guess you have to be
% selecting the neighbors? i forget if i made this into indices



