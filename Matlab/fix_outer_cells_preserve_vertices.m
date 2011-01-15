% Erodes the image by deleting pixels close to the boundary
% from all 4 sides. Does not delete a pixel if it is vertex
% The code here is horrible. it could be much
% faster and probably otherwise better. But that's ok, people don't use
% this function much anyway ;)


function bords = fix_outer_cells_preserve_vertices(x)


Ys = size(x, 1);
Xs = size(x, 2);

% the old definition of vertices, different from bwmorph('branchpoints')
f = @(x) (x(2,2) && (sum(x(:)) >= 4));
lut = makelut(f, 3);
vertplot = applylut(x, lut);



bords = x;

for i = 1:Ys
    for j = 1:Xs
        if x(i, j)
            if ~vertplot(i, j)
                bords(i, j) = 0;
            end
            break;
        end
    end
end

for i = 1:Ys
    for j = Xs:-1:1
        if x(i, j)
            if ~vertplot(i, j)
                bords(i, j) = 0;
            end
            break;
        end
    end
end




for i = 1:Xs
    for j = 1:Ys
        if x(j, i)
            if ~vertplot(j, i)
                bords(j, i) = 0;
            end
            break;
        end
    end
end

for i = 1:Xs
    for j = Ys:-1:1
        if x(j, i) 
            if ~vertplot(j, i)
                bords(j, i) = 0;
            end
            break;
        end
    end
end