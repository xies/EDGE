% given a string like "22-100", return the numbers 22 and 100
% also accepts string like "555" in which case second is returned as null

function [out] = parse_range(in)

out = [NaN NaN];

if sum(in=='-') == 0
    out(1) = str2double(in);
    out(2) = str2double(in);
    return;
end


nums = regexp(in, '[0-9]+');
if length(nums) ~= 2
    return
end
n1 = nums(1);
n2 = nums(2);

out(1) = str2double( in(n1:n2-2) );
out(2) = str2double( in(n2:end)  ); 