function out = my_mean(in)

% mean that ignores NaN values
% only operates along the 2nd dimension

out = zeros(size(in,1),1);

for i = 1:size(in, 1)
    x=in(i,:);
    x=x(~isnan(x));
    out(i) = mean(x);
end
    