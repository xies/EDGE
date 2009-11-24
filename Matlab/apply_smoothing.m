function out = apply_smoothing(data, handles, NAME)

dt = handles.info.seconds_per_frame / 60;
sigma = get(handles.smoothing_strength_slider, 'Value');

[smo.single, hf.single, roc.single] = ...
    smooth_rochange_2(data.single, sigma,dt);   

if iscell(data.neighbors_1st)
    for i = 1:length(data.neighbors_1st)
        [smo.neighbors_1st{i}, hf.neighbors_1st{i}, roc.neighbors_1st{i}] = ...
        smooth_rochange_2(data.neighbors_1st{i}, sigma,dt);  
    end
else
    [smo.neighbors_1st, hf.neighbors_1st, roc.neighbors_1st] = ...
        smooth_rochange_2(data.neighbors_1st, sigma,dt);   
end
if iscell(data.neighbors_2nd)
    for i = 1:length(data.neighbors_2nd)
        [smo.neighbors_2nd{i}, hf.neighbors_2nd{i}, roc.neighbors_2nd{i}] = ...
            smooth_rochange_2(data.neighbors_2nd{i}, sigma,dt);   
    end
else
    [smo.neighbors_2nd, hf.neighbors_2nd, roc.neighbors_2nd] = ...
        smooth_rochange_2(data.neighbors_2nd, sigma,dt);   
end

if get(handles.button_smoothed, 'Value')
    out.single = smo.single;
    out.neighbors_1st = smo.neighbors_1st;
    out.neighbors_2nd = smo.neighbors_2nd;
elseif get(handles.button_high_frequency, 'Value')
    out.single = hf.single;
    out.neighbors_1st = hf.neighbors_1st;
    out.neighbors_2nd = hf.neighbors_2nd;
elseif get(handles.button_rate_of_change, 'Value')
    if strcmp(NAME, 'Area (borders)') ...
            || strcmp(NAME, 'Area (polygon)') ...
            || strcmp(NAME, 'Area (orthogonal)')
        % constriction rate
        out.single = -roc.single;
        out.neighbors_1st = -roc.neighbors_1st;
        out.neighbors_2nd = -roc.neighbors_2nd;
    else
        out.single = roc.single;
        out.neighbors_1st = roc.neighbors_1st;
        out.neighbors_2nd = roc.neighbors_2nd;
    end
end