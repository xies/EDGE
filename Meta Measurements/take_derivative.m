function out = take_derivative(in, windowsize)

out = diff(smooth(in, windowsize)).';