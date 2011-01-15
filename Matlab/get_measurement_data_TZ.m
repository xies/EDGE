% Find the values of the selected measurement from time T0 to T1 and from 
% depth Z0 to Z1 (gets the measurement directly from the downdown menu)

function [data name unit] = get_measurement_data_TZ(...
    handles, dropdown_val, cell_index)%, T0, T1, Z0, Z1)

drop_str = handles.allmeasures;
big_measure_name = drop_str{dropdown_val};
% big_measure_name(big_measure_name == ':') = '-';

% then, get the right part of it 
data = handles.loaded_measurements.(genvarname(big_measure_name)).data(:, :, cell_index);
name = handles.loaded_measurements.(genvarname(big_measure_name)).name;
unit = handles.loaded_measurements.(genvarname(big_measure_name)).unit;
