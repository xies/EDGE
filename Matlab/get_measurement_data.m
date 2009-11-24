function [data name units] = get_measurement_data(handles, channel, ...
    measure_filename, T, Z, storedT, storedZ, C)
% this is the function that reads from stored_properties and finds the
% correct measurement. 

% INPUT ARGUMENTS:
% channel           -- the name of the channel, e.g., 'Membranes'
% measure_filename  -- e.g., 'ellipse_properties'
% index             -- the index of which measurement you want within that measure_filename, e.g., 1
% [T Z]             -- the T and Z coordinate we want
% [storedT storedZ] -- the T and Z but translated for the indices of the stored_properties struct array, i.e., starting from 1 and increasing
% C                 -- the desired cell index

% OUTPUT ARGUMENTS:
% data              -- the data itself, either a scalar or an array
% name              -- the name of that measurement, a string
% units             -- the units of that measurement, a string


% built_in_measurements = {'Area', 'Perimeter', 'Centroid-x', 'Centroid-y'};
if isempty(channel)
% if it's one of the built-in measurments
    name = measure_filename;
    cell = handles.embryo.getCellGraph(T, Z).getCell(C);
    switch measure_filename
        case 'Area'
            if isempty(cell)
                data = NaN;
            else
                data = cell.area * handles.info.microns_per_pixel^2;
            end
            units= 'micron^2';
        case 'Perimeter'
            if isempty(cell)
                data = NaN;
            else
                data = cell.perimeter * handles.info.microns_per_pixel;
            end
            units= 'microns';
        case 'Centroid-x'
            if isempty(cell)
                data = NaN;
            else
                data = cell.centroid * handles.info.microns_per_pixel;
                data = data(2);
            end
            units= 'microns';
        case 'Centroid-y'
            if isempty(cell)
                data = NaN;
            else
                data = cell.centroid * handles.info.microns_per_pixel;
                data = data(1);
            end
            units= 'microns';
        otherwise
            disp('Error in get_measurement_data.m: measurement name does not contain :: identifier but is also not a built-in measurement');
            data = [];
            name = [];
            units= [];
            return;
    end
else

    % extract the channel and measurement name from 'measure'
%     IDENTIFIER = '::';
%     dots = strfind(measure, IDENTIFIER);
%     channel = measure(1:dots(1)-1);
%     filename = measure(dots(1)+length(IDENTIFIER):dots(2)-1);
%     measurename = measure(dots(2)+length(IDENTIFIER):end);

    data  = handles.stored_properties.(channel).(measure_filename).data{storedT, storedZ, C};
    name  = handles.stored_properties.(channel).(measure_filename).names;
    units = handles.stored_properties.(channel).(measure_filename).units;

end

