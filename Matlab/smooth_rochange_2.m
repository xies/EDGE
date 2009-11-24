function [smo_out,hf_out,roc_out]=smooth_rochange_2(data, sigma, dt)
% smooths each column of the incoming data, 
% using a Gaussian 
% smoothing with window size 5 and sigma input.

% returns...
% the smoothed data (smo_out)
% the high frequency component, original - smoothed, (hf_out)
% the rate of changes based on the smoother data (roc_out)

if isempty(data)
    smo_out = [];
    hf_out = [];
    roc_out = [];
    return;
end

if sigma > 0
    WINDOW_SIZE = 10;%1 + floor(1*sigma);

    % convert sigma from seconds to index units
    stdev = sigma / dt; 
    rowdata = data.';

    % pad the data
    left_pad = repmat(rowdata(:, 1), 1, WINDOW_SIZE);
    right_pad = repmat(rowdata(:, end), 1, WINDOW_SIZE);
    background = [left_pad rowdata right_pad];

    % smooth using the built-in "smoothts" function with the Gaussian option
    smoothed = smoothts(background, 'g', WINDOW_SIZE, stdev);

    % unpad (extract it from the background elements)
    smo_out = smoothed(:, 1+WINDOW_SIZE:end-WINDOW_SIZE);
    smo_out = smo_out.'; % back to original orientation (columns)
else
    smo_out = data;
end

% calculate high frequency
hf_out = smo_out - data;

% find the rate of change
% (diff acts along 1st dimention by default)
roc_out = diff(smo_out) / dt;
% add NaN to the end so it is the same length
roc_out(end+1, :) = NaN;

