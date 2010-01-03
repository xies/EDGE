function import_copymove_files(handles, image_file, dest, rm_files, sample_image)
% copy the imported files from the import directory to the EDGE directory
% image_file is the image file *function* for this data but with the source
    % already specified so that only a time and layer are needed
% dest is the destination directory
% rm_files is boolean and says whether or not to remove the originals
% sample_image is just a sample image to check if it needs rgb2gray


for time_i = handles.info.start_time:handles.info.end_time
    set(handles.text_processing_time,  'String', num2str(time_i));
    for layer_i = handles.info.bottom_layer:my_sign(handles.info.top_layer-handles.info.bottom_layer):handles.info.top_layer
        if get(handles.radiobutton_stop, 'Value')
            set(handles.radiobutton_stop, 'Value', 0);
            readyproc(handles, 'ready');
            msgbox('Import interrupted by user - please try again', 'Import failed');    
            guidata(hObject, handles);
            return;
        end

        set(handles.text_processing_layer, 'String', num2str(layer_i));
        drawnow

        copymove_src = image_file(time_i, layer_i);
        if rm_files
            movefile(copymove_src, dest);
        else
            copyfile(copymove_src, dest);
        end
    end
end

% make sure all the images are grayscale

if ndims(sample_image) == 3
    readyproc(handles, 'finalizing');
    for time_i = handles.info.start_time:handles.info.end_time
        set(handles.text_processing_time,  'String', num2str(time_i));
        for layer_i = handles.info.bottom_layer:my_sign(handles.info.top_layer-handles.info.bottom_layer):handles.info.top_layer
            set(handles.text_processing_layer, 'String', num2str(layer_i));
            drawnow

            filename = handles.info.image_file(time_i, layer_i, dest);
            out = imread(filename);
            in = rgb2gray(out);
            imwrite(in, filename, handles.file_ext);
         end
    end
end