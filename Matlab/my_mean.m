function out = my_mean(in)
% takes the mean along the 2nd dimension and ignores NaN values

out = zeros(size(in, 1), 1);
for i = 1:size(in, 1)
    row = in(i, :);
    row(isnan(row)) = [];  % get rid of NaN values
    out(i) = mean(row);
end

