function [img r] = drawCellSmall(embryo, t, z, c)
% draws an image of the cell in a bounding box. the image is returned in
% "img" and the coordinates of the top-left corner are returned in "r"

vcoords = round(embryo.getCellGraph(t, z).getCell(c).vertexCoords);
r = min(vcoords);

vcoords(:, 1) = vcoords(:, 1) - min(vcoords(:, 1));
vcoords(:, 2) = vcoords(:, 2) - min(vcoords(:, 2));
img = poly2mask(vcoords(:, 2), vcoords(:, 1), ceil(max(vcoords(:,1))), ceil(max(vcoords(:, 2))));
