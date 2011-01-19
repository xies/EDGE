function out = my_sum(in)

% sum that ignores NaN values
% only operates along the 2nd dimension

out = zeros(size(in,1),1);

for i = 1:size(in, 1)
    x=in(i,:);
    x=x(~isnan(x));
    out(i) = sum(x);
end
    