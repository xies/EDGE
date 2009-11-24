function out = imagestripe(Ys, Xs)

out = zeros(Ys, Xs);
out(floor(Ys/2)-1:ceil(Ys/2)+1, :) = 1;
