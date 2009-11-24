function handle = write_image_filename_function(SRC, name, z_posn, z_digits, t_posn, t_digits, fixed, data_set_name)

fid = fopen(fullfile(SRC, 'image_filename.m'), 'w');

fprintf(fid, '%s\n', 'function data = image_filename(time_i, layer_i, src)');
fprintf(fid, '%s\n', '% ** This is an automatically generated function');
fprintf(fid, '%s%s%s\n', '% ** created at ', date_and_time, ' by write_image_filename_function.m');
fprintf(fid, '%s\n%s\n\n', '% ** Inputs the time, layer, and source directory of a data set.', ...
    '% ** Outputs the filename of that image.');
fprintf(fid, '%s%s\n\n', '% ** For data set: ', data_set_name);
fprintf(fid, '%s%s%s\n\n', 'filename = ''', name, ''';');
fprintf(fid, '%s%u%s\n', 'z_name = sprintf(strcat(''%.'', num2str(', z_digits, ...
    '), ''u''), layer_i);');
fprintf(fid, '%s%u%s%u%s%u%s\n\n', 'filename(', z_posn, ':', z_posn, '+', z_digits, ...
    '-1) = z_name;');
if ~fixed
    fprintf(fid, '%s%u%s\n', 't_name = sprintf(strcat(''%.'', num2str(', t_digits, ...
        '), ''u''), time_i);');
    fprintf(fid, '%s%u%s%u%s%u%s\n\n', 'filename(', t_posn, ':', t_posn, '+', t_digits, ...
        '- 1) = t_name;');
end
fprintf(fid, '%s\n', 'data = fullfile(src, filename);');
fclose(fid);

current_dir = pwd;
cd(SRC);
handle = @image_filename;
cd(current_dir);

% ** this function should make a function that looks like this: **

% function data = image_filename(time_i, layer_i, src)
% filename = 'Series001_t000_z0_ch00.tif';
% 
% z_name = sprintf(strcat('%.', num2str(3), 'u'), layer_i);
% filename(8 : 8 + 3 - 1) = z_name;
% 
% t_name = sprintf(strcat('%.', num2str(2), 'u'), time_i);
% filename(4 : 4 + 2 - 1) = t_name;
% 
% data = fullfile(src, filename);


