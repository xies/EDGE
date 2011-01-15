% Find the vertices given a border image "bords". Returns them in an array.

function [vertlist] = find_vertices(bords)


% old implementation (slower): 
% f = @(x) (x(2,2) && (sum(x(:)) >= 4));
% lut = makelut(f, 3);
% vertplot = applylut(bords, lut);

vertplot = bwmorph(bords, 'branchpoints');
%{
but the above is different- it gives "pre-merged" vertices
actually, this might be much better. i should think about this,
might save time or whatever. or if i don't merge vertices... then all the recursive
stuff would have just been a waste of time (yup...) well.....
well-- sometimes it's useful if you don't want small edges, like for 
measuring neighbor changes?
%}


[y x] = find(vertplot);
vertlist = [y(:) x(:)];