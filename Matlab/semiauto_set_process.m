function semiauto_set_process(handles)

% if in bandpass mode, turn them visible
if get(handles.radiobutton_process_bandpass, 'Value')
    
    set(handles.text_bandpass_high,'Visible','on');
    set(handles.text_bandpass_low,'Visible','on');
    set(handles.info_text_bandpass_high,'Visible','on');
    set(handles.info_text_bandpass_low,'Visible','on');
    set(handles.info_text_preprocessing_threshold,'Visible','on');
    set(handles.text_preprocess_thresh,'Visible','on');
    
    set(handles.info_text_disperse_threshold,'Visible','off');
    set(handles.text_disperse_threshold,'Visible','off');
    
else
    
    set(handles.text_bandpass_high,'Visible','off');
    set(handles.text_bandpass_low,'Visible','off');
    set(handles.info_text_bandpass_high,'Visible','off');
    set(handles.info_text_bandpass_low,'Visible','off');
    set(handles.info_text_preprocessing_threshold,'Visible','off');
    set(handles.text_preprocess_thresh,'Visible','off');
    
    set(handles.info_text_disperse_threshold,'Visible','on');
    set(handles.text_disperse_threshold,'Visible','on');
    
end

end

% this is screwy-- should only add cells by adding edges
% it doesn't work because if i cut into a cell to add a new one, the
% vertices of that existing cell don't get updated. basically, this
% function does not take into account the effects it has on existing cells,
% but just throws a new one in the mix. AddEdge does this properly with no
% problems
% set(handles.vec_add_cell, 'Enable', 'off');