function [location whichbutton] = get_mouse_location_zoom(handles)

loc = get(handles.axes1, 'CurrentPoint');
whichbutton =get(handles.figure1, 'SelectionType');


location(1) = loc(1,2);
location(2) = loc(1,1);

