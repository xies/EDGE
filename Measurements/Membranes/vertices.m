function [data names units] = vertices(embryo, getMembranes, t, z, c, dx, dz, dt, other)
% computes the ellipse properties at (t, z) for Cell i given the Embryo4D,
% the membranes, and the relevant resolutions

% the names
names{1} = 'Vertex-x';
names{2} = 'Vertex-y';
names{3} = '# of vertices';
names{4} = '# of neighbors';

% the units
units{1} = 'microns';
units{2} = 'microns';
units{3} = '';
units{4} = '';

verts = embryo.getCell(c, t, z).vertexCoords;

data{1} = verts(:, 2) * dx;
data{2} = verts(:, 1) * dx;
data{3} = length(verts);
data{4} = length(embryo.getCell(c, t, z).neighbors);
