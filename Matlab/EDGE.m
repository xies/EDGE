function varargout = EDGE(varargin)
    % EDGE M-file for EDGE.fig
    %      EDGE, by itself, creates a new EDGE or raises the existing
    %      singleton*.
    %
    %      H = EDGE returns the handle to a new EDGE or the handle to
    %      the existing singleton*.
    %
    %      EDGE('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in EDGE.M with the given input arguments.
    %
    %      EDGE('Property','Value',...) creates a new EDGE or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before EDGE_OpeningFunction gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to EDGE_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help EDGE

    % Last Modified by GUIDE v2.5 19-Dec-2009 23:23:33

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @EDGE_OpeningFcn, ...
                       'gui_OutputFcn',  @EDGE_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT
 

% --- Executes just before EDGE is made visible.
function EDGE_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to EDGE (see VARARGIN)
    clc;
    
    % info from semiauto, if you got here by switching
    if ~isempty(varargin) && length(varargin) > 1
        passinfo = varargin{2};
    else
        passinfo = [];
    end
    
    % if we are starting in the wrong directory, fix it...
    main_dir = mfilename('fullpath');
    this_filename = mfilename;
    main_dir = main_dir(1:end-length(this_filename));
    cd(main_dir);
    
    %%% read the names of the Data Sets!!!!
    % but only include those that have been exported in the list
    datanames = get_folder_names(fullfile('..', 'DATA_GUI'));
    embryonames = strcat('..', filesep, 'DATA_GUI', filesep, datanames, filesep, 'embryo_data.mat');
    gooddatanames = cell(0);
    for i = 1:length(embryonames)
        if exist(embryonames{i}, 'file')
            gooddatanames{length(gooddatanames)+1} = datanames{i};
        end
    end
    
    if ~isempty(gooddatanames)
        %%% set the dropdown string
        set(handles.dropdown_datasets, 'String', gooddatanames);
    

        % if you press a key, run keyFunction (defined by me)
        have_KeyPressFcn = findobj('KeyPressFcn', '');
        for i = 1:length(have_KeyPressFcn)
            set(have_KeyPressFcn(i), 'KeyPressFcn', @keyFunction);
        end

        % default data set
        if ~isempty(passinfo)
            default_data_set = passinfo.data_set;
            % set the dropdown to be at the right thing
            val = find(strcmp(gooddatanames, passinfo.data_set));
            set(handles.dropdown_datasets, 'Value', val);    
            % make sure we're in the right directory
            cd(main_dir);
        else
            default_data_set_number = 1;
            default_data_set = gooddatanames{default_data_set_number};  % just take the first one
            set(handles.dropdown_datasets, 'Value', default_data_set_number);
        end
        
        
        
        
        handles = clear_data_set_maingui(handles, default_data_set);    
        
        
        % if there are arguments passed in (i.e., if we switched from EDGE)
        % then set them here
        if ~isempty(passinfo)
            set(handles.t_slider, 'Value', passinfo.t_slider);
            t_slider_Callback(hObject, [], handles);
            set(handles.z_slider, 'Value', passinfo.z_slider);
            z_slider_Callback(hObject, [], handles);
        end

    else
        msgbox(strcat('There are currently no data sets which have been exported by semiauto.', ...
            'To continue, first process a data set with semiauto'), 'EDGE: no data sets', 'error');
    end
    
    % initialize the image (after the data box is checked!)
    set(hObject,'toolbar','figure');

%     set(handles.figure1,'CloseRequestFcn',@closeGUI)  
    set(handles.figure1,'WindowButtonUpFcn',@mouseFunction);
%     set(handles.axes1, 'ButtonDownFcn', @mouseFunction);
%     set(handles.axes1, 'NextPlot', 'add');
    
%     set(handles.axes1, 'Units', 'pixels');
%     set(handles.figure1, 'Units', 'pixels');
%     set(handles.main_panel, 'Units', 'pixels');
    set(handles.axes2, 'Color', 'none');
    
    % Choose default command line output for EDGE
    handles.output = hObject;

    set(hObject,'toolbar','figure');
    % Update handles structure
 
%     set(handles.neighbors_panel,'SelectionChangeFcn',@neighbor_panel_change);
    set(handles.smoothing_panel,'SelectionChangeFcn',@smoothing_panel_change);
    set(handles.panel_3d_z_or_t, 'SelectionChangeFcn', @panel_3d_change);
    set(handles.panel_neighbors_averages, 'SelectionChangeFc', @averages_panel_change);
    set(handles.radiobutton_panel_3d,'SelectionChangeFcn', @wireframe_surface_panel_change);
%     set(handles.export_panel,        'SelectionChangeFcn',@export_panel_selection);
    
    guidata(hObject, handles);
    
    % UIWAIT makes EDGE wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
 

% --- Outputs from this function are returned to the command line.
function varargout = EDGE_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;
 

% --- Executes on slider movement.
function z_slider_Callback(hObject, eventdata, handles)
    % hObject    handle to z_slider (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'Value') returns position of slider
    %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    %obtains the slider value from the slider component
    sliderValue = get(handles.z_slider,'Value');
    sliderValue = round(sliderValue);
    set(handles.z_slider, 'Value', sliderValue);
    
    %puts the slider value into the edit text component
    showValue = handles.info.bottom_layer - sliderValue * ...
        sign(handles.info.bottom_layer - handles.info.top_layer);

    set(handles.z_text,'String', num2str(showValue));

    neighbor_val = str2double(get(handles.edit_neighbors_order, 'String'));
    if neighbor_val > 0
        handles = maingui_initialize_neighbors(handles, neighbor_val);
    end
    
    handles = slider_callbacks_draw_image_slice(handles);
    slider_callbacks_draw_3D_cell(handles);
    slider_callbacks_draw_measurement(handles);
    
%     handles = slider_callbacks(handles, [1 1 1]);
    
%     set(hObject,'toolbar','figure');
    % Update handles structure
    guidata(hObject, handles);
 

% --- Executes during object creation, after setting all properties.
function z_slider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
 


% --- Executes on slider movement.
function t_slider_Callback(hObject, eventdata, handles)
    % hObject    handle to t_slider (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'Value') returns position of slider
    %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    sliderValue = get(handles.t_slider,'Value');
    sliderValue = round(sliderValue);
    set(handles.t_slider, 'Value', sliderValue);
    
    showValue = handles.info.start_time + sliderValue * ...
        sign(handles.info.end_time - handles.info.start_time);
    
    %puts the slider value into the edit text component
    set(handles.t_text,'String', num2str(showValue));

    neighbor_val = str2double(get(handles.edit_neighbors_order, 'String'));
    if neighbor_val > 0
        handles = maingui_initialize_neighbors(handles, neighbor_val);
    end
    
    handles = slider_callbacks_draw_image_slice(handles);
    slider_callbacks_draw_3D_cell(handles);
    slider_callbacks_draw_measurement(handles);
    
    % Update handles structure
    guidata(hObject, handles);
 

% --- Executes during object creation, after setting all properties.
function t_slider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
 


% --- Executes on button press in cbox_raw.
function cbox_raw_Callback(hObject, eventdata, handles)
    val1 = get(handles.cbox_raw, 'Value');
    val2 = get(handles.cbox_bord,'Value');
    val3 = get(handles.cbox_poly,'Value');
    val4 = get(handles.cbox_other, 'Value');
    
    if ~(val1 || val2 || val3 || val4)
        set(handles.cbox_raw,'Value', 1);
    else
    handles = slider_callbacks_draw_image_slice(handles);
    
%         handles = slider_callbacks(handles, [1 0 0]);
        guidata(hObject, handles);
    end
 

% --- Executes on button press in cbox_bord.
function cbox_bord_Callback(hObject, eventdata, handles)
    val1 = get(handles.cbox_raw, 'Value');
    val2 = get(handles.cbox_bord,'Value');
    val3 = get(handles.cbox_poly,'Value');
    val4 = get(handles.cbox_other, 'Value');
    
    if ~(val1 || val2 || val3 || val4)
        set(handles.cbox_bord,'Value', 1);
    else
%         handles = slider_callbacks(handles, [1 0 0]);
        handles = slider_callbacks_draw_image_slice(handles);
        guidata(hObject, handles);
    end
 

% --- Executes on button press in cbox_poly.
function cbox_poly_Callback(hObject, eventdata, handles)
    val1 = get(handles.cbox_raw, 'Value');
    val2 = get(handles.cbox_bord,'Value');
    val3 = get(handles.cbox_poly,'Value');
    val4 = get(handles.cbox_other, 'Value');
    
    if ~(val1 || val2 || val3 || val4)
        set(handles.cbox_poly,'Value', 1);
    end
    
    if ~get(handles.cbox_poly,'Value')
        handles = clear_cell(handles);
        handles = slider_callbacks_draw_image_slice(handles);
        slider_callbacks_draw_3D_cell(handles);
        slider_callbacks_draw_measurement(handles);
    else
        handles = slider_callbacks_draw_image_slice(handles);
    end

    guidata(hObject, handles);
    
    
function cbox_other_Callback(hObject, eventdata, handles)
    val1 = get(handles.cbox_raw,  'Value');
    val2 = get(handles.cbox_bord, 'Value');
    val3 = get(handles.cbox_poly, 'Value');
    val4 = get(handles.cbox_other,'Value');
    
    if ~(val1 || val2 || val3 || val4)
        set(handles.cbox_other, 'Value', 1);
        return;
    end
    
    % if you're turning it on
    if val4
        % if there's only one to choose from, automatically select that
        if length(handles.channelnames) == 1
            handles.activeChannels = 1;
        else
            % choose the set of channels
            [selection ok] = listdlg('ListString', handles.channelnames, 'Name', 'Select channel(s)');
            if ~ok
                return
            end
            handles.activeChannels = selection;
        end
    else % if you're turning it off
        handles.activeChannels = [];
    end

    handles = slider_callbacks_draw_image_slice(handles);
    
    guidata(hObject, handles);


function keyFunction(src,evnt)
    %keyPressFcn automatically takes in two inputs
    %src is the object that was active when the keypress occurred
    %evnt stores the data for the key pressed

    %brings in the handles structure in to the function
    handles = guidata(src);

    k = evnt.Key; %k is the key that is pressed

%     disp(k);

    switch k
        case 'semicolon'  % time slider to the left by 1
    %         pause(0.01) %allows time to update
            %define hObject as the object of the callback that we are going to use
            %in this case, we are mapping the enter key to the add_pushbutton
            %therefore, we define hObject as the add pushbutton
            %this is done mostly as an error precaution
            
            
            hObject = handles.t_slider; 

            sliderValue = get(handles.t_slider,'Value');
            if sliderValue-1 >= get(handles.t_slider, 'Min')
                set(handles.t_slider, 'Value', sliderValue-1);

                %call the add pushbutton callback.  
                %the middle argument is not used for this callback
                t_slider_Callback(hObject, [], handles);
            end
        case 'quote'  % time slider to the right by 1
            hObject = handles.t_slider;

            sliderValue = get(handles.t_slider,'Value');
            if sliderValue+1 <= get(handles.t_slider, 'Max')
                set(handles.t_slider, 'Value', sliderValue+1);
                t_slider_Callback(hObject, [], handles);
            end
        case 'leftbracket'   % z  slider to the left by 1
            hObject = handles.z_slider;

            sliderValue = get(handles.z_slider,'Value');
            if sliderValue-1 >= get(handles.z_slider, 'Min')        
                set(handles.z_slider, 'Value', sliderValue-1);
                z_slider_Callback(hObject, [], handles);
            end
        case 'rightbracket'  % z slider to the right by 1
            hObject = handles.z_slider;

            sliderValue = get(handles.z_slider,'Value');

            if sliderValue+1 <= get(handles.z_slider, 'Max')        
                set(handles.z_slider, 'Value', sliderValue+1);
                z_slider_Callback(hObject, [], handles);
            end
        case 'q'  %'escape'  % quit the gui
            closeGUI = handles.figure1;
            close(closeGUI);
        otherwise
          % nothing
    end
    
% for the data cursor stuff - failure    
% function output_txt = dcupdate(obj, eventdata)
%     disp('hi')
%     output_txt = '';
%     pos = get(eventdata,'Position')
        
function mouseFunction(hObject,evnt)
    handles = guidata(hObject);    
    
    % only select a cell if the polygon mode is on
    polyv = get(handles.cbox_poly, 'Value');
    if ~polyv
        return
    end
     
    Ys = handles.info.Ys; Xs = handles.info.Xs;
    location = get_mouse_location_zoom(handles);
    
    % make sure we're in the figure
    if location(1) > Ys || location(1) < 1 || location(2) > Xs || location(2) < 1
        return
    end
    
      
    [T Z] = getTZ(handles);
    
    newCell = handles.embryo.getCellGraph(T, Z).activeCellAtPoint(location);
    if isempty(newCell) % clicked outside the polygons
        return;
    end
    
    if get(handles.button_manually_select_cells, 'Value')
    % if we are selecting manually such that we can increase the size of
    % handles.activeCell indefinitely. note that in this case
    % handles.activeCell is a cell array
        addCell = 1;
        % if you click on the same cell it unselects
        if ~isempty(handles.activeCell)
            for i = 1:length(handles.activeCell)
              if handles.activeCell(i) == newCell.index;
                  handles.activeCell(i) = [];
                  addCell = 0;
                  break;
              end
            end
        end
        if addCell
            handles.activeCell(length(handles.activeCell)+1) = newCell.index;
        end
        
        
    else
    % in the other case where only one cell at a time can be selected, such
    % that handles.activeCell is either equal to a Cell object or []
        if handles.activeCell == newCell.index % unselect if it's the same one
            handles = clear_cell(handles);
        else
            handles.activeCell = newCell.index;
            
            % if there are neighbors, add them on
            val = str2double(get(handles.edit_neighbors_order, 'String'));
            if val > 0  
                handles = maingui_initialize_neighbors(handles, val);
            end
            
        end
        
    end
    

        
    % deal with the drawing stuff. if you've selected more that one
    % cell, then the celltext box needs to show '-' or something like
    % that, and also it won't draw in 3d and you can't select neighbors
    if isempty(handles.activeCell)
        handles = clear_cell(handles);
    elseif ~get(handles.button_manually_select_cells, 'Value')  % single cell selection
        set(handles.cell_text, 'String', num2str(handles.activeCell)); 
        slider_callbacks_draw_measurement(handles);
        slider_callbacks_draw_3D_cell(handles);
    else % multiple cell selection
%         handles = clear_cell(handles);
        slider_callbacks_draw_measurement(handles);
    end
        
    % redraw
%     handles = slider_callbacks_draw_image_slice_dots_maingui(handles);
    handles = slider_callbacks_draw_image_slice(handles);
    
    guidata(hObject, handles);
 



function smoothing_panel_change(hObject, eventdata)
    %retrieve GUI data, i.e. the handles structure
    handles = guidata(hObject); 
%     handles = slider_callbacks(handles, [0 0 1]);
    slider_callbacks_draw_measurement(handles);
        
    %updates the handles structure
    guidata(hObject, handles);
 
    
function averages_panel_change(hObject, eventdata)
    handles = guidata(hObject); 
    slider_callbacks_draw_measurement(handles);
        
    %updates the handles structure
    guidata(hObject, handles);
    
function panel_3d_change(hObject, eventdata)
    handles = guidata(hObject);
    % onyl allow to take pictures for temporal movies
    if get(handles.radiobutton_3d_temporal, 'Value')
        set(handles.button_make_movie, 'String', 'take picture');
        set(handles.movie_fps_text, 'Visible', 'off');
        set(handles.movie_fps, 'Visible', 'off')
    else
        set(handles.button_make_movie, 'String', 'make movie');
        set(handles.movie_fps_text, 'Visible', 'on');
        set(handles.movie_fps, 'Visible', 'on')
    end
    
    slider_callbacks_draw_3D_cell(handles);
    guidata(hObject, handles);

% --- Executes on selection change in dropdown_measurements.
function dropdown_measurements_Callback(hObject, eventdata, handles)
    % if needed, load the file
    dropdown_val = get(handles.dropdown_measurements, 'Value');
    big_measure_name = handles.allmeasures{dropdown_val};
    
    already_loaded = fieldnames(handles.loaded_measurements);
    if ~any(strcmp(already_loaded, genvarname(big_measure_name))) 
set(handles.text_readyproc, 'Visible', 'on');
set(handles.text_readyproc, 'String', 'Loading...');
set(handles.text_readyproc, 'ForegroundColor', [1 0 0]);
drawnow;
        load(fullfile(handles.src.measurements, [big_measure_name '.mat']));  % loads 'data', 'name', 'unit'
        handles.loaded_measurements.(genvarname(big_measure_name)).data = data;
        handles.loaded_measurements.(genvarname(big_measure_name)).name = name;
        handles.loaded_measurements.(genvarname(big_measure_name)).unit = unit;
set(handles.text_readyproc, 'Visible', 'off');
    end
    guidata(hObject, handles);


    slider_callbacks_draw_measurement(handles);
 

% --- Executes during object creation, after setting all properties.
function dropdown_measurements_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
 

% --- Executes on selection change in dropdown_datasets.
function dropdown_datasets_Callback(hObject, eventdata, handles)
    all_datasets = get(handles.dropdown_datasets, 'String');
    new_dataset_number = get(handles.dropdown_datasets, 'Value');
    new_dataset = all_datasets{new_dataset_number};

    if strcmp(new_dataset, handles.data_set)
        return;  % if you selected the same one, do nothing
    end

    handles = clear_data_set_maingui(handles, new_dataset);
    guidata(hObject, handles);
    
 

% --- Executes during object creation, after setting all properties.
function dropdown_datasets_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
 

% --- Executes on slider movement.
function smoothing_strength_slider_Callback(hObject, eventdata, handles)
    sliderValue = get(handles.smoothing_strength_slider,'Value');

%     sliderValue = round(sliderValue);
%     set(handles.smoothing_strength_slider,'Value', sliderValue);
    
    %puts the slider value into the edit text component
    set(handles.smoothing_strength_text,'String', num2str(sliderValue));

    slider_callbacks_draw_measurement(handles);
%     handles = slider_callbacks(handles, [0 0 1]);
    guidata(hObject, handles);
 

% --- Executes during object creation, after setting all properties.
function smoothing_strength_slider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end   
 

% function closeGUI(hObject, eventdata)
%     handles = guidata(hObject);
%     clear handles.activeCell;
%     clear handles.embryo;
%     clear handles;
%     delete(gcf);
 

% --- Executes on button press in export_button.
function export_button_Callback(hObject, eventdata, handles)
    if isempty(handles.activeCell)
        msgbox('Error: No cells selected.', 'Export aborted', 'error');
%         uiwait;
        return;
    end

    % create the directories if necessary
    dirname = fullfile(handles.program_dir, 'DATA_OUTPUT');
    [a b c] = mkdir(dirname);
    dirname = fullfile(dirname, handles.data_set);
    [a b c] = mkdir(dirname);
%     dirname = fullfile(dirname, date_and_time_saving);
%     [a b c] = mkdir(dirname);

    [T Z] = getTZ(handles);
    cell_indices = handles.activeCell;
    save(fullfile(dirname, date_and_time_saving), 'cell_indices', 'T', 'Z');
    
    msgbox('Successfully exported data to the DATA_OUTPUT folder.');
    
%     % select which data sets to export
%     datanames = get_folder_names(fullfile('..', 'DATA_GUI'));
%     default_val = find(strcmp(datanames, handles.data_set));
%     if length(datanames) > 1
%         [selection_data ok] = listdlg('ListString', datanames, 'Name', 'Select data sets for export', ...
%             'ListSize', [300 300], 'InitialValue', default_val);
%         if ~ok
%             return;
%         end
%     else
%         selection_data = default_val;
%     end
%     data_set = handles.data_set;
%     selectbutton_val = get(handles.button_manually_select_cells, 'Value');
%     
%     if length(selection_data) == 1
%         % note: if you don't select any cells, it just exports for all
%         % cells..!
%         
%         % change datasets immediately so that it gives the correct list of
%         % measurements
%         if ~strcmp(data_set, datanames{selection_data})
%             handles = clear_data_set_maingui(handles, datanames{selection_data});
%             handles.activeCell = handles.embryo.getCellGraph(handles.info.master_time, handles.info.master_layer).activeCellIndices;
%             set(handles.button_manually_select_cells, 'Value', 1);
%         end
%             
%                 
%         % for just one data set, do the usual thing
%         default_val = get(handles.dropdown_measurements, 'Value');
%         [selection_meas ok] = listdlg('ListString', handles.allmeasures, ...
%             'Name', 'Select measurements sets for processing', ...
%             'InitialValue', default_val, 'ListSize', [300 300]);
%         if ~ok
%             return;
%         end
%                 
%         
%     else
%         res = questdlg('For exporting multiple data sets, all cells and all measurements will be exported.', ...
%                   'Export all cells all measurements?', 'Ok', 'Cancel', 'Ok');
%         if strcmp(res, 'Cancel')
%             return;
%         end
%    
%    
% %         % for multiple data sets, pick the measurements
% %         [measurementchannelsall measurementnamesall] = get_measurement_file_names_specific_ch(...
% %             fullfile(handles.program_dir, 'Measurements'), ...
% %             handles.all_channelnames);
% %         allmeasurements_list = strcat(measurementchannelsall, '::', measurementnamesall);
% %         allmeasurements_list = [handles.builtin; allmeasurements_list(:)];
% %         % let the user slect
% %         [selection_meas_multi ok] = listdlg('ListString', allmeasurements_list, ...
% %             'Name', 'Select measurements sets for export', 'ListSize', [300 300]);
% %         if ~ok
% %             return;
% %         end
% 
%     end
%     
%     
%     
%      % old:
%     % first, ask which measurements to export (by FILE?? no, no point, you can select multiple)
% %     [measurementchannelsall measurementnamesall] = get_measurement_file_names(handles);
% %     built_in_measurements = {'Area', 'Perimeter', 'Centroid-x', 'Centroid-y'};
% % 
% %     listnames = strcat(measurementchannelsall, '::', measurementnamesall);
% %     listnames = [built_in_measurements listnames];
% 
%     
% 
%     for export_datasets_i = 1:length(selection_data)
%         % if it's already at this data set, don't need to clear_data_set
%         if ~strcmp(handles.data_set, datanames{selection_data(export_datasets_i)})
%             handles = clear_data_set_maingui(handles, datanames{selection_data(export_datasets_i)});
%         end
%         
%         % if no cells selected, select all of them
%         if isempty(handles.activeCell)
%             handles.activeCell = handles.embryo.getCellGraph(handles.info.master_time, handles.info.master_layer).activeCellIndices;
%             set(handles.button_manually_select_cells, 'Value', 1);
%         end
%         
%         set(handles.text_readyproc, 'Visible', 'on');
%         set(handles.text_readyproc, 'String', 'Exporting...');
%         set(handles.text_readyproc, 'ForegroundColor', [1 0 0]);
%         drawnow;
%     
%         
%         if length(selection_data) > 1
%             selection_meas = 1:length(handles.allmeasures);
%         end
%         
% %         % choose the appropriate measurements if multiple data sets
% %         if length(selection_data) > 1
% %             % now, find all the measurements that are common to 
% %             % measurementnames_thisdata AND measurementnames_multi
% %             % the former is a list of all measurements for this data set,
% %             % and the latter is the selection of all choices
% %             selection_meas = CStrAinBP(handles.allmeasures, allmeasurements_list);
% %         end
%         
%         
%     
%         % create the directories if necessary
%         dirname = fullfile(handles.program_dir, 'DATA_OUTPUT');
%         [a b c] = mkdir(dirname);
%         dirname = fullfile(dirname, handles.data_set);
%         [a b c] = mkdir(dirname);
%         dirname = fullfile(dirname, date_and_time_saving);
%         [a b c] = mkdir(dirname);
% 
%         x_vals = [];
%         savedata = [];
%         savedata_neighbors = [];
% 
%         % for each measurement you selected
%         for export_measurements_i = 1:length(selection_meas)
% 
%             % can make use of the file get_measurement_data_TZ
%             [T Z] = getTZ(handles);
% 
%             % plots the data
%     %         totaldata = [];  % can define size of this later
% 
%             % get the measure name from the dropdown menu
%             dropdown_val = selection_meas(export_measurements_i);
% 
% 
%             if handles.fixed
%                 x_vals = handles.info.bottom_layer:my_sign(handles.info.top_layer-handles.info.bottom_layer):handles.info.top_layer;
%                 x_vals = x_vals * handles.info.microns_per_z_step;
%                 x_vals = x_vals(:);
% 
%                 if get(handles.button_manually_select_cells, 'Value')  % multiple cells selected
% 
%                     totaldata = cell(abs(handles.info.top_layer-handles.info.bottom_layer)+1, length(handles.activeCell));
%                     for i = 1:length(handles.activeCell)
%                         [data name units] = get_measurement_data_TZ(...
%                             handles, dropdown_val, handles.activeCell(i), handles.info.start_time, handles.info.end_time, ...
%                             handles.info.bottom_layer, handles.info.top_layer);           
%     %                     data = convert_cell_data_to_numerical(data);
%                         totaldata(:, i) = data(:);
%                     end
% 
%                     savedata = totaldata; 
%     %                 if get(handles.radiobutton_neighbors_averages, 'Value') % if we need to average
%     %                     savedata = my_mean(savedata);
%     %                 end
% 
%                 else  % single cell selected
%                     [data name units] = get_measurement_data_TZ(...
%                         handles, dropdown_val, handles.activeCell, handles.info.start_time, handles.info.end_time, ...
%                         handles.info.bottom_layer, handles.info.top_layer);   
%     %                 data = convert_cell_data_to_numerical(data);
%                     savedata = data;
% 
%                     % get the measurements for the neighbors as well
%                     if ~isempty(handles.activeCellNeighbors)    
%                         savedata_neighbors = cell(length(handles.activeCellNeighbors), 1);
%                         for i = 1:length(handles.activeCellNeighbors)
%                             for j = 1:length(handles.activeCellNeighbors{i})
%                                 [data name units] = get_measurement_data_TZ(...
%                                     handles, dropdown_val, handles.activeCellNeighbors{i}(j), handles.info.start_time, handles.info.end_time, ...
%                                     handles.info.bottom_layer, handles.info.top_layer);           
%     %                             data = convert_cell_data_to_numerical(data);
%                                 totaldata(:, j) = data(:);
%                             end
%                             if get(handles.radiobutton_neighbors_averages, 'Value') % if we need to average
%                                 savedata_neighbors{i} = my_mean(totaldata);
%                             else
%                                 savedata_neighbors{i} = totaldata;
%                             end
%                         end
%                     end
% 
%                 end
% 
%             else  % time series
% 
%                 x_vals = handles.info.start_time:handles.info.end_time;
%                 x_vals = x_vals * handles.info.seconds_per_frame / 60;
%                 x_vals = x_vals(:);
% 
%                 if get(handles.button_manually_select_cells, 'Value')  % for multiple cells
% 
%                     totaldata = cell(abs(handles.info.end_time-handles.info.start_time)+1, ...
%                                    abs(handles.info.top_layer-handles.info.bottom_layer)+1, ...
%                                    length(handles.activeCell));
%                     for i = 1:length(handles.activeCell)
% 
%                         [data name units] = get_measurement_data_TZ(...
%                             handles, dropdown_val, handles.activeCell(i), handles.info.start_time, handles.info.end_time, ...
%                             handles.info.bottom_layer, handles.info.top_layer);
%     %                     data = convert_cell_data_to_numerical(data);
%                         totaldata(:, :, i) = data;
%                     end
% 
%     %                 if get(handles.radiobutton_neighbors_averages, 'Value') % if we need to average
%     %                     savedata = my_mean(totaldata);
%     %                 else
%                         savedata = totaldata;
%     %                 end
%                 else  % just one cell selected
% 
%                     % in this case, it's a bit tricky. if there are no neighbors, we do
%                     % something special and plot all layers, with red-->blue
%                     % representing top-->bottom. however, if there are neighbors we
%                     % skip this and just plot the middle cell at this current depth,
%                     % and the same goes for all the neighbors
%                     if isempty(handles.activeCellNeighbors)
%                         [data name units] = get_measurement_data_TZ(...
%                             handles, dropdown_val, handles.activeCell, handles.info.start_time, handles.info.end_time, ...
%                             handles.info.bottom_layer, handles.info.top_layer);
%     %                     data = convert_cell_data_to_numerical(data);
% 
%                         % if we are on averages, we average over all layers!!
%     %                     if get(handles.radiobutton_neighbors_averages, 'Value')
%     %                         savedata = my_mean(data);
%     %                     else
%                             savedata = data;
%     %                     end
%                     else
%                         % in the neighbors case, we first need to plot the active
%                         % (middle) cell in the non-special way (i.e., we just want this
%                         % layer, not all layers), and then we plot all the neighbors as
%                         % well.
%                         [data name units] = get_measurement_data_TZ(...
%                             handles, dropdown_val, handles.activeCell, handles.info.start_time, handles.info.end_time, ...
%                             Z, Z);   
%     %                     data = convert_cell_data_to_numerical(data);
%                         savedata = data;
% 
%                         % now we plot for all the neighbors
%                         savedata_neighbors = cell(length(handles.activeCellNeighbors), 1);
%                         for i = 1:length(handles.activeCellNeighbors)
%                             totaldata = cell(abs(handles.info.end_time-handles.info.start_time)+1, length(handles.activeCellNeighbors{i}));
%                             for j = 1:length(handles.activeCellNeighbors{i})
%                                 [data name units] = get_measurement_data_TZ(...
%                                     handles, dropdown_val, handles.activeCellNeighbors{i}(j), handles.info.start_time, handles.info.end_time, ...
%                                     Z, Z);           
%     %                             data = convert_cell_data_to_numerical(data);
%                                 totaldata(:, j) = data(:);
%                             end
%     %                         if get(handles.radiobutton_neighbors_averages, 'Value') % if we need to average
%     %                             savedata_neighbors{i} = my_mean(totaldata);
%     %                         else
%                                 savedata_neighbors{i} = totaldata;
%     %                         end
%                         end
%                     end
% 
% 
%                 end
% 
%             end
% 
%             % saving
%             name = handles.allmeasures{selection_meas(export_measurements_i)};
%             name(name == ':') = '-';
% 
%     %         IDENTIFIER = '::';
%     %         dots = strfind(name, IDENTIFIER);
%     %         dots = dots(end);
%     %         if ~isempty(dots)
%     %             name = name(dots + length(IDENTIFIER) : end);
%     %         end
% 
%             % save cell numbers
%             cell_indices = handles.activeCell;
% 
%             filename = fullfile(dirname, name);
%             save(filename, 'x_vals', 'savedata', 'savedata_neighbors', 'cell_indices');
% 
%         end
%     end  % end for multiple data sets
%     if ~strcmp(handles.data_set, data_set)
%         handles = clear_data_set_maingui(handles, data_set);
%     end
%     set(handles.button_manually_select_cells, 'Value', selectbutton_val);
% 
%     set(handles.text_readyproc, 'Visible', 'off');
%     msgbox(strcat('Successfully exported data to the DATA_OUTPUT folder.'));
    


function cell_text_Callback(hObject, eventdata, handles)
    ind = str2double(get(handles.cell_text, 'String'));
    if ind > 0 && round(ind)==ind
        [T Z] = getTZ(handles);
        handles.activeCell = ind;
        slider_callbacks_draw_measurement(handles);
        slider_callbacks_draw_3D_cell(handles);
%         handles = slider_callbacks_draw_image_slice_dots_maingui(handles);
        handles = slider_callbacks_draw_image_slice(handles);
        guidata(hObject, handles);
    else
        handles.activeCell = [];
        handles = clear_cell(handles);
        msgbox('Cell number must be a positive integer', '', 'error');
    end
    guidata(hObject, handles);


function cell_text_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end



function axis_min_Callback(hObject, eventdata, handles)
    min = str2double(get(handles.axis_min, 'String'));
    max = str2double(get(handles.axis_max, 'String'));
    if min >= max
        set(handles.axis_min, 'String', '-');
        msgbox('min must be less than max', '', 'error');
    elseif isnan(min)
        set(handles.axis_min, 'String' ,'-');
%         handles = slider_callbacks(handles, [0 0 1]);
        slider_callbacks_draw_measurement(handles);
    else
%         handles = slider_callbacks(handles, [0 0 1]);
        slider_callbacks_draw_measurement(handles);
    end
%     guidata(hObject, handles);   

function axis_min_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


function axis_max_Callback(hObject, eventdata, handles)
    min = str2double(get(handles.axis_min, 'String'));
    max = str2double(get(handles.axis_max, 'String'));
    if min >= max
        set(handles.axis_min, 'String', '-');
        msgbox('max must be greater than min', '', 'error');
    elseif isnan(min)
        set(handles.axis_min, 'String' ,'-');
%         handles = slider_callbacks(handles, [0 0 1]);
        slider_callbacks_draw_measurement(handles);
    else
%         handles = slider_callbacks(handles, [0 0 1]);
        slider_callbacks_draw_measurement(handles);
    end
%     guidata(hObject, handles);

function axis_max_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


    
function wireframe_surface_panel_change(hObject, eventdata)
    %retrieve GUI data, i.e. the handles structure
    handles = guidata(hObject);   
    if handles.activeCell > 0
%         handles = slider_callbacks(handles, [0 1 0]);
        slider_callbacks_draw_3D_cell(handles);
    end
    guidata(hObject, handles);


% function export_panel_selection(hObject, eventdata)
%     handles = guidata(hObject);   
%     switch get(eventdata.NewValue,'Tag')
%         case 'radiobutton_export_all_layers'
%             handles.export_all_layers = 1;
%         case 'radiobutton_export_this_layer'
%             handles.export_all_layers = 0;
%         otherwise
%             %
%     end
%     guidata(hObject, handles);


function button_make_movie_Callback(hObject, eventdata, handles)
    % create the directories if necessary
    dirname = fullfile(handles.program_dir, 'DATA_MOVIES');
    [a b c] = mkdir(dirname);
    dirname = fullfile(dirname, handles.data_set);
    [a b c] = mkdir(dirname);
    
    % shoud we make a movie (or a picture?)
    movie_on = ~(handles.fixed || get(handles.radiobutton_3d_temporal, 'Value'));
    
    filename = fullfile(dirname, strcat('Cell', num2str(handles.activeCell), ...
        '_Created_', date_and_time_saving));
    
    if movie_on
        % create the AVI file
        fps = str2double(get(handles.movie_fps, 'String'));
        aviobj = avifile(filename, 'fps', fps);
    end
 
 
% this code is useless- can just uses the axis
    axesprops = get(handles.axes2);
    names = fieldnames(axesprops);
 
    axes(handles.axes2);
    axis_props = axis;
% the axis properties of the current image. will use this for all images. 
    f = figure(99);
    ax = gca;
    set(ax, 'Color', 'none');
%     axes(ax);
    
       
%     indTime = 0;
    for time_i = handles.info.start_time:handles.info.end_time
        cla;
        % set z = NaN so it doesn't highlight any layer
        if get(handles.radiobutton_3d_spatial, 'Value')
            slider_callbacks_draw_3D_cell(handles, time_i, NaN); 
        else
            [T Z] = getTZ(handles);
            slider_callbacks_draw_3D_cell(handles, NaN, Z); 
        end
        
        axis(axis_props);
        
        for i = 1:length(names)
            cont = strfind(names{i}, 'Camera');
            if ~isempty(cont) && cont(1) == 1
                set(ax, names{i}, axesprops.(names{i}));
            end
        end

        if movie_on
            F = getframe(f);
            aviobj = addframe(aviobj,F);
        end
       
%         indTime = indTime + 1;
    end

    if ~movie_on
        saveas(f, filename, 'tif');
    else
        aviobj = close(aviobj);
    end
    close(f);
    
    
   
function movie_fps_Callback(hObject, eventdata, handles)


function movie_fps_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cbox_show_membs_Callback(hObject, eventdata, handles)
%     handles = slider_callbacks(handles, [0 1 0]);
    slider_callbacks_draw_3D_cell(handles);
%     guidata(hObject, handles);

function cbox_show_surface_Callback(hObject, eventdata, handles)
%     handles = slider_callbacks(handles, [0 1 0]);
    slider_callbacks_draw_3D_cell(handles);
%     guidata(hObject, handles);

function cbox_show_vertices_Callback(hObject, eventdata, handles)
%     handles = slider_callbacks(handles, [0 1 0]);
%     guidata(hObject, handles);
    slider_callbacks_draw_3D_cell(handles);

function cbox_show_centroids_Callback(hObject, eventdata, handles)
%     handles = slider_callbacks(handles, [0 1 0]);
%     guidata(hObject, handles);    
    slider_callbacks_draw_3D_cell(handles);

function cbox_show_centroid_fit_Callback(hObject, eventdata, handles)
%     handles = slider_callbacks(handles, [0 1 0]);
%     guidata(hObject, handles);
    slider_callbacks_draw_3D_cell(handles);

function cbox_show_other_channels_Callback(hObject, eventdata, handles)
%     handles = slider_callbacks(handles, [0 1 0]);
%     guidata(hObject, handles);
    slider_callbacks_draw_3D_cell(handles);

function button_switch_to_semiauto_Callback(hObject, eventdata, handles)
    handles = guidata(hObject);
    
    passinfo.data_set = handles.data_set;
    passinfo.t_slider = get(handles.t_slider, 'Value');
    passinfo.z_slider = get(handles.z_slider, 'Value');
    
    % extract all the camera properties
    axesinfo = [];
%     axesprops = get(handles.axes1);
%     names = fieldnames(axesprops);
%     for i = 1:length(names)
%         cont = strfind(names{i}, 'Camera');
%         if ~isempty(cont) && cont(1) == 1
%             axesinfo.(names{i}) = axesprops.(names{i});
%         end
%     end
    
    delete(handles.figure1);
    semiauto('dummy', passinfo, axesinfo);


function handles = edit_neighbors_order_Callback(hObject, eventdata, handles)
    val = str2double(get(handles.edit_neighbors_order, 'String'));
    if isempty(val) || isnan(val) || val < 0 || round(val) ~= val
        set(handles_edit_neighbors_order, 'String', '0');
%         set(handles.radiobutton_3d_spatial, 'Enable', 'off');
%         set(handles.radiobutton_3d_temporal, 'Enable', 'off');
        return;
    end
    
    
    % if a cell is selected at the time, then replot. otherwise no need to
    if ~isempty(handles.activeCell)
        handles = maingui_initialize_neighbors(handles, val);
%         handles = slider_callbacks_draw_image_slice_dots_maingui(handles); 
        handles = slider_callbacks_draw_image_slice(handles);
        slider_callbacks_draw_measurement(handles);    
        if get(handles.radiobutton_3d_temporal, 'Value')
            slider_callbacks_draw_3D_cell(handles);
        end
    end
    
    
    guidata(hObject, handles);
    
    
    


function edit_neighbors_order_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function button_manually_select_cells_Callback(hObject, eventdata, handles)
    if get(handles.button_manually_select_cells, 'Value')
        handles.activeCell = [];
        set(handles.button_manually_select_cells, 'String', 'unselect all');
        set(handles.edit_neighbors_order, 'Enable', 'off');
        set(handles.button_neighbors_plus, 'Enable', 'off');
        set(handles.button_neighbors_minus, 'Enable', 'off');
    else
        handles.activeCell = [];
        set(handles.button_manually_select_cells, 'String', 'select manually');
        set(handles.edit_neighbors_order, 'Enable', 'on');
        set(handles.button_neighbors_plus, 'Enable', 'on');
        set(handles.button_neighbors_minus, 'Enable', 'on');
    end
%     handles = clear_cell(handles);
    handles = slider_callbacks_draw_image_slice(handles);
    slider_callbacks_draw_3D_cell(handles);
    slider_callbacks_draw_measurement(handles);    
    guidata(hObject, handles);


function button_select_all_Callback(hObject, eventdata, handles)
    % activate the manual mode
    set(handles.button_manually_select_cells, 'Value', 1);
    set(handles.button_manually_select_cells, 'String', 'unselect all');
    set(handles.edit_neighbors_order, 'Enable', 'off');
    set(handles.button_neighbors_plus, 'Enable', 'off');
    set(handles.button_neighbors_minus, 'Enable', 'off');
        
    handles = clear_cell(handles);
    [T Z] = getTZ(handles);
    handles.activeCell = handles.embryo.getCellGraph(handles.info.master_time, handles.info.master_layer).activeCellIndices;
%     handles = slider_callbacks_draw_image_slice_dots_maingui(handles); 
    handles = slider_callbacks_draw_image_slice(handles);
    slider_callbacks_draw_measurement(handles);
    guidata(hObject, handles);
   


function button_neighbors_plus_Callback(hObject, eventdata, handles)
    val = str2double(get(handles.edit_neighbors_order, 'String'));
    set(handles.edit_neighbors_order, 'String', num2str(val+1));
    handles = edit_neighbors_order_Callback(hObject, eventdata, handles);
    guidata(hObject, handles);
    
function button_neighbors_minus_Callback(hObject, eventdata, handles)
    val = str2double(get(handles.edit_neighbors_order, 'String'));
    if val == 0
        return;
    end
    set(handles.edit_neighbors_order, 'String', num2str(val-1));
    handles = edit_neighbors_order_Callback(hObject, eventdata, handles);
    guidata(hObject, handles);
