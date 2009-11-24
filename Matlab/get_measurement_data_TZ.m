function [data name units] = get_measurement_data_TZ(...
    handles, dropdown_val, cell_index, T0, T1, Z0, Z1)
% find the values of the selected measure from time T0 to T1 and from depth Z0
% to Z1 (gets the measurement directly from the downdown menu)

% % get the measure name from the dropdown menu
% dropdown_val = get(handles.dropdown_measurements, 'Value');
% drop_str = get(handles.dropdown_measurements, 'String');
drop_str = handles.allmeasures;
big_measure_name = drop_str{dropdown_val};

% split it into channel, filename, and measurement name
IDENTIFIER = '::';
dots = strfind(big_measure_name, IDENTIFIER);
if isempty(dots)  % for one of the built-in measurements
    measure_channel = [];
    measure_filename = big_measure_name;
    % if it's one of the built in measurements, set the channel to empty
    % and store the name in measure_filename
else
    measure_channel = big_measure_name(1:dots(1)-1);
    measure_filename = big_measure_name(dots(1)+length(IDENTIFIER):dots(2)-1);
%     measure_specific_name = big_measure_name(dots(2)+length(IDENTIFIER):end);
    % if this measure file has several measurements associated with it, we need
    % to find out which one we want. in other words, data, name, and units
    % might be arrays, and we need to know which element to pick out. to do
    % this, we call get_measurement_names and see while one the measurename
    % matches up with
%     allmeasurementnames = get_measurement_names(handles, channel,
%     measure_filename);
    allindices = strfind(drop_str, measure_filename);
    % get the first non-empty value -- this is where all the measurements
    % for that measure_filename start
    for i = 1:length(allindices)
        if ~isempty(allindices{i})
            start_ind = i;
            break;
        end
    end
    index = find(strcmp(drop_str, big_measure_name));
    index = index - start_ind + 1;
end





% % built-in measurement
% if isempty(measure_channel)
    data = cell(abs(T1-T0)+1, abs(Z1-Z0)+1);
    % get the data at all x_vals and depths
    indTime = 1;
    for time_i = T0:my_sign(T1-T0):T1
        indLayer = 1;
        for layer_i = Z0:my_sign(Z1-Z0):Z1
            [data_out] = ...
                get_measurement_data(handles, measure_channel, measure_filename, ...
                time_i, layer_i, abs(time_i - handles.info.start_time) + 1, ...
                abs(layer_i - handles.info.bottom_layer) + 1, cell_index); 
            if ~isempty(measure_channel)
                % get the right one if there are several measurements
                % this does not apply to built-in measurements

                if iscell(data_out)
                    data{indTime, indLayer} = data_out{index};
                else
                    data{indTime, indLayer} = data_out(index);
                end

            else
                data{indTime, indLayer} = data_out;
            end

            indLayer = indLayer + 1;
        end
        indTime = indTime + 1;
    end




% else
%     indT0 = abs(T0 - handles.info.start_time) + 1;
%     indT1 = abs(T1 - handles.info.start_time) + 1;
%     indZ0 = abs(Z0 - handles.info.bottom_layer) + 1;
%     indZ1 = abs(Z1 - handles.info.bottom_layer) + 1;
%     data  = handles.stored_properties.(measure_channel).(measure_filename).data(...
%         indT0:my_sign(indT1-indT0):indT1, ...
%         indZ0:my_sign(indZ1-indZ0):indZ1, ...
%         cell_index);
%     name  = handles.stored_properties.(measure_channel).(measure_filename).names;
%     units = handles.stored_properties.(measure_channel).(measure_filename).units;
% 
%     name  = name{index};
%     units = units{index};
% 
%     
% end






% do it this way so that we get the names/units for the built-in props as well
[dummy name units] = get_measurement_data(handles, measure_channel, measure_filename, ...
    handles.info.master_time, handles.info.master_layer, ...
    abs(handles.info.master_time-handles.info.start_time)+1, ...
    abs(handles.info.master_layer-handles.info.bottom_layer)+1, 1);
if ~isempty(measure_channel)
%     name = handles.stored_properties.(measure_channel).(measure_filename).names;
%     units = handles.stored_properties.(measure_channel).(measure_filename).names;
    name  = name{index};
    units = units{index};
end
% actually we already know the correct name from the selection in the
% dropdown menu, but we do this to get the units
