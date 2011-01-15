% Controls the text in semiauto that says whether the data set is exported.

function exported_text_controller(handles, val)

switch val
    case 'exported'
        set(handles.text_exported, 'String', 'Exported');
        set(handles.text_exported, 'ForegroundColor', [0 1 0]);
        set(handles.button_switch_to_explorer, 'Enable', 'on');
    
    case 'not exported'
        set(handles.text_exported, 'String', 'Not exported');
        set(handles.text_exported, 'ForegroundColor', [1 1 0]);
        set(handles.button_switch_to_explorer, 'Enable', 'off');
end

drawnow;
    