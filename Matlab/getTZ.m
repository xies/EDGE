% Extract the current T and Z from handles.

function [T Z] = getTZ(handles)
T = str2double(get(handles.t_text,'String'));
Z = str2double(get(handles.z_text,'String'));

% t = get(handles.t_slider,'Value');
% z = get(handles.z_slider,'Value');

if isnan(handles.info.seconds_per_frame) % fixed data set
%     t = 0;
    T = 0;
end

% t and z are the indeces starting from 0
% and are used as indeces to access the java arrays

% T and Z are the corresponding values for the numbers of the
% images themselves, and are used with read_file

% note: this means that, in this version, the true slider "value"
% and what is being shown are different. this is because a slider
% must always move UP as it moves RIGHT, but i don't always want this.
% luckily the user only seeing the text box, so it works well like this.