% Asks if you want to save changes if there are any changes.
% This function is to be called whenever you are about to leave a data set, 
% either by changing to another one or closing the GUI. 


function handles = exit_data_set_semiauto(handles)

fields = changed_data_info(handles);
% ask about unsaved changes to DATA_INFO
if ~isempty(fields)  % if there are changes
    othermsg = 'There are unsaved changes to DATA_INFO in the following fields:';
    msg2 = 'Save them now?';
    res = questdlg([othermsg fields msg2], 'Saving changes', ...
                 'Save changes', 'Discard changes', 'Save changes');  
    if strcmp(res, 'Save changes')
        write_data_info(handles);
%         othermsg = 'Successfully updated the following fields in DATA_INFO.csv:';
%         msg = [othermsg fields];
%         msgbox(msg, 'Save successful');
%         save_embryo(handles);
    elseif strcmp(res, 'Discard changes')
        % we need to save the embryo to be consistent with the old stuff,
        % so load it back in and save over
        handles.info = read_data_info(handles.data_set);
        handles = update_embryo(handles);
    end
end


% ask about unsaved changes to the embryo data
if handles.embryo.changed
    msg = 'There are unsaved changes to the Embryo data.';
    res = questdlg(msg, 'Saving changes', ...
                 'Save changes', 'Discard changes', 'Save changes');  
    if strcmp(res, 'Save changes')
        % this is the only time we call save_embryo (!!)
        readyproc(handles, 'saving');
        handles = save_embryo(handles);
    elseif strcmp(res, 'Discard changes')
        % do nothing
    end
end




