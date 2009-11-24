function out = convert_cell_data_to_numerical(data)
% the user-defined measurements return data as a cell array like 
% data{1} = ... , data{2} = ... ,etc. the function
% get_measurement_data_TZ finds the correct one of these and returns a cell
% array with dimensions t() x z(), where each element is the correct data
% point. but we don't know if the data itself is a single numerical value
% or an array. for plottng, we just want to get the single numerical value
% (say, the first). but for exporting we want all of them. since we plan to
% use the get_measurement_data_TZ function for exporting (it's admittedly
% wasteful but deals with built-in stuff)


out = zeros(size(data));
for i = 1:size(data, 1)
    for j = 1:size(data, 2)
        if isnumeric(data{i, j}(1))
            out(i, j) = data{i, j}(1);
        else 
            out(i, j) = NaN;
        end
    end
end
    
    