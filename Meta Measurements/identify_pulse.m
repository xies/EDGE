clear all; close all; clc;

% measurement that identifies apical constriction pulses

data_set = 'Histone032';
cell_inds = []; % all cells
measurement = 'Membranes--basic_2d--Area';
layers_from_top = 3;

[data_all data_apical] = extract_measurement(data_set, measurement, cell_inds, layers_from_top);

% looking at the apical data only
% first, compute the rate of change of area
rate = diff(data_apical, ????);

% then, identify pulses by thresholding this rate at 1 st above meant
pulses = rate > mean(rate) + std(rate);







