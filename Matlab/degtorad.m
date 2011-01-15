% Convert degrees to radians.

function rad = degtorad(deg)
% because the MathWorks is weird and made this function part of some
% toolbox. also, they changed the name which was annoying.

rad = (pi/180)*deg;