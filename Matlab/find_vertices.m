function [vertlist] = find_vertices(bords)
% just gets the simple vertices by definition and returns the list

f = @(x) (x(2,2) && (sum(x(:)) >= 4));
lut = makelut(f, 3);

vertplot = applylut(bords, lut);

[y x] = find(vertplot);

vertlist = [y(:) x(:)];



%{
could just do

verplot = bwmorph(bords, 'branchpoints');


but this is different- it gives "pre-merged" vertices
actually, this might be much better. i should think about this,
might save time or whatever. or if i don't merge vertices... then all the recursive
stuff would have just been a waste of time..



%}