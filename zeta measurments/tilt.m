function [data names units] = tilt(embryo, getMembranes, t, z, c, dx, dz, dt)
% computes the tilt angles at (t, z) for Cell i given the Embryo4D,
% the membranes, and the relevant resolutions
% based on linear regresion in the vicinity of z

% the names
names{1} = 'Tilt-x';
names{2} = 'Tilt-y';
names{3} = 'Tilt-tot';
names{4} = 'Tilt-dir';

% the units
units{1} = 'degree';
units{2} = 'degree';
units{3} = 'degree';
units{4} = 'degree';

% paramters used in this routine
n_sel = round(6/dz);   % number of z layers chosen to calculate tilt


% extract coordinates of centroids
centroids = Cell.centroidStack(embryo.getCellStack(c, t)) * dx;
y_values = centroids(:, 1);
x_values = centroids(:, 2);
z_values = (1:length(x_values)) * dz;
z_values = z_values(:);



% select coordinates of 'n_sel' closest layers to z that are not NaN
ind_an = find(~isnan(x_values));
z_dist = abs(z_values(embryo.translateZ(z)+1)-z_values(ind_an)); 
[q,ind] = sort(z_dist);
ind = ind(1:min(n_sel,length(ind)));   % use 'n_sel' layers or less 
x_sel = x_values(ind_an(ind));
y_sel = y_values(ind_an(ind));
z_sel = z_values(ind_an(ind));


% compute tilt via principle components
X = [x_sel,y_sel,z_sel];
[coeff,score,roots] = princomp(X);    
 
% first component, vector pointing towards largest variation of centroids
v = coeff(:,1);  
if (v(3)<0)
    v=-v;
end

ang_x = atand(v(1)/v(3));   % angles in x
ang_y = atand(v(2)/v(3));   % and y

% spherical coordinates
ang_tot = atand(sqrt(v(1)^2+v(2)^2)/abs(v(3))); % angle of tilt with z axis
ang_dir = atan2(v(2),v(1))/pi*180;              % and angle in xy plane

% output

data{1} = ang_x;
data{2} = ang_y;
data{3} = ang_tot;
data{4} = ang_dir;

