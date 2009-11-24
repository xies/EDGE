function handles = semiauto_info_text_callback(handles, name)
button_name = strcat('info_text_', name);

string = get(handles.(button_name), 'String');
var = str2double(string);

old = my_num2str(handles.info.(name));
if strcmp(string, '-') || strcmp(string, 'NaN')
    set(handles.(button_name), 'String', '-');
    handles.info.(name) = NaN;
elseif isnan(var)
    set(handles.(button_name), 'String', old);  
else
    handles.info.(name) = var;
end