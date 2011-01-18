function drop = find_nucleus_drop(position)

drop_percent = 0.8;

x = position(:);

if isnan(my_mean(x(1:3)'))
    drop = NaN;
else
    drop = find(x < my_mean(x(1:3)')*drop_percent, 1);
end

if isempty(drop)
    drop = NaN;
end

