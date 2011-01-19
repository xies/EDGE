function out = take_derivative(in)

out = diff(smooth(in, 3)).';