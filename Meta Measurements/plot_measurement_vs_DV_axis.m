function plot_measurement_vs_DV_axis(data_set, measurement, cell_inds, name)


[data_allz data_topz] = extract_measurement(data_set, measurement, cell_inds);


%%
figure;
plot(my_mean(data_allz));
xlabel('D-V');
ylabel(measurement)