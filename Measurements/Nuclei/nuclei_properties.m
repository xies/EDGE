function [data names units] = nuclei_properties(embryo, getNuclei, t, z, c, dx, dz, dt,others)

% the names
names{1} = 'intensity';
names{2} = 'position-x';
names{3} = 'position-y';
names{4} = 'position-z';
names{5} = 'lower end';
names{6} = 'upper end'; 
names{7} = 'length';
names{8} = 'depth';
% position x, y , z, /absolute relative distance to top

% the units
units{1} = 'intensity';
units{2} = 'micron';
units{3} = 'micron';
units{4} = 'micron';
units{5} = 'micron';
units{6} = 'micron';
units{7} = 'micron';
units{8} = 'micron';

data = num2cell(NaN(length(units), 1));


%% optional: process only given cell

% if c ~= 353;% 209 %353
% % if (c <= 358 || c>369)
%      return;
% end



%% do for every layer
nuclei = getNuclei(t, z);
[cell_img R] = drawCellSmall(embryo, t, z, c); % draw the Cell

nuclei = nuclei(R(1):R(1)+size(cell_img, 1)-1, R(2):R(2)+size(cell_img, 2)-1);
nuclear_intensity = sum(sum(nuclei(cell_img)));

data{1} = nuclear_intensity;


%% whole-cell properties

% get layers that were tracked; all others are NaN in x_values 
centroids = Cell.centroidStack(embryo.getCellStack(c, t)) * dx;
x_values = centroids(:, 2);




if (z ~= embryo.masterLayer ||  sum(~isnan(x_values))<20)
     return;
end

c



%% in case we are in the master layer
nuclear_intensity_profile = zeros(embryo.z, 1);
ind_all = ~isnan(x_values);   % tracked and repaired cells
nuclei = cell(embryo.z, 4);
% loop over all layers to get...
for z_i = embryo.lowestTracked(c, t) : my_sign(embryo.highestTracked(c, t) - embryo.lowestTracked(c, t)) : embryo.highestTracked(c, t)
    if (~isnan(x_values(embryo.translateZ(z_i)+1)))
        % ...intensity profile for fitting
        [cell_img R] = drawCellSmall(embryo, t, z_i, c); % draw the Cell
        nucleus = getNuclei(t, z_i);
        nucleus = nucleus(R(1):R(1)+size(cell_img, 1)-1, R(2):R(2)+size(cell_img, 2)-1);
        if (z_i==embryo.lowestTracked(c, t))
           avg_bg = mean(mean(nucleus(cell_img)));
           nuclear_intensity = 0;
        else
            %nucleus = nucleus>3*avg_bg;
        end
        nuclear_intensity = mean(mean(nucleus(cell_img)));
        nuclear_intensity_profile(embryo.translateZ(z_i)+1) = nuclear_intensity;
        % ...all values of nuclei labelling and coordinates for center of mass
        nuclei{embryo.translateZ(z_i)+1,1} = double(nucleus);
        nuclei{embryo.translateZ(z_i)+1,2} = cell_img;
        [X,Y] = meshgrid(R(2):R(2)+size(cell_img, 2)-1, R(1):R(1)+size(cell_img, 1)-1);
        nuclei{embryo.translateZ(z_i)+1,3} = X*dx;
        nuclei{embryo.translateZ(z_i)+1,4} = Y*dx;
    else
        % if not tracked, take previous nucleus
        nuclei{embryo.translateZ(z_i)+1,1} = double(nucleus);
        nuclei{embryo.translateZ(z_i)+1,2} = cell_img;
        nuclei{embryo.translateZ(z_i)+1,3} = X*dx;
        nuclei{embryo.translateZ(z_i)+1,4} = Y*dx;
        ind_all(z_i) = 1;
    end
    
end


%% Fitting sum of two Gaussians to intensity profile to extract position and width of nucleus

% z-axis for fitting
z_values = (1:length(x_values)) * dz;
z_values = z_values(:);
ind = find(~isnan(x_values));
% initialvalues of fit parameters
[max_val max_ind] = max(nuclear_intensity_profile);
max_z = z_values(max_ind);
std_z = std(z_values(ind));
% set fitting options
fitopt = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[max_val/1.5,max_z-10,std_z/3, max_val/5/1.5,max_z-10,std_z*10/3],...
    'Upper',[max_val*1.5,max_z+10,std_z*3, max_val/5*1.5,max_z+10,std_z*10],...
    'Startpoint',[max_val,max_z,std_z, max_val/5,max_z,std_z*10]);
% perform fit, fitfuntion is sum of two Gaussians
cfun = fit(z_values(ind),nuclear_intensity_profile(ind),'gauss2',fitopt);

%     figure(1);clf;
%     plot(z_values,nuclear_intensity_profile(1:end),'o',...
%          z_values,feval(cfun,z_values),'-r');
% max_ind = embryo.unTranslateZ(max_ind);


%% calculate center of mass of nucleus (x and y position)

% estimate range over which to look for nucleus
ind_nuc = find(z_values>=cfun.b1-cfun.c1*2 & z_values<=cfun.b1+cfun.c1*2 & ...
    ind_all==1);
norm_nuc = zeros(1,1);
com_nuc = zeros(2,1);
for i=ind_nuc(1):ind_nuc(end)
    if (i==107)
    end
    nucleus = nuclei{i,1};
    cell_img = nuclei{i,2};
    X = nuclei{i,3};
    Y = nuclei{i,4};
    norm_nuc = norm_nuc + sum(sum(nucleus(cell_img)));
    com_nuc(1) = com_nuc(1) + sum(sum(double(nucleus(cell_img)).*X(cell_img)));
    com_nuc(2) = com_nuc(2) + sum(sum(double(nucleus(cell_img)).*Y(cell_img)));
end
com_nuc = com_nuc/norm_nuc;

%% absolute and relative nucleus length and distance from cell apex


%centroid coordinats
x_val = x_values;
y_val = centroids(:, 1);
z_val = z_values;

% remove NaN values
for i = length(x_val):-1:1
    if isnan(x_val(i))
        x_val(i) = [];
        y_val(i) = [];
        z_val(i) = [];
    end
end

% cell length
cell_length = sum(sqrt(diff(x_val).^2 + diff(y_val).^2 + diff(z_val).^2));

% depth of nucleus
ind_top = find(z_val > cfun.b1);
x_v = [com_nuc(1);x_val(ind_top)];
y_v = [com_nuc(2);y_val(ind_top)];
z_v = [cfun.b1;z_val(ind_top)];
nucleus_depth = sum(sqrt(diff(x_v).^2 + diff(y_v).^2 + diff(z_v).^2));

% length of nucleus
ind_nuc = find(z_val>cfun.b1-cfun.c1*2 & z_val<cfun.b1+cfun.c1*2);
x_v = [x_val(ind_nuc(1));x_val(ind_nuc);x_val(ind_nuc(end))];
y_v = [y_val(ind_nuc(1));y_val(ind_nuc);y_val(ind_nuc(end))];
z_v = [cfun.b1-cfun.c1*2;z_val(ind_nuc);cfun.b1+cfun.c1*2];
nucleus_length = sum(sqrt(diff(x_v).^2 + diff(y_v).^2 + diff(z_v).^2));
%% output
data{2} = com_nuc(1);   % position of nucleaus in x
data{3} = com_nuc(2);   % position of nucleaus in y
data{4} = cfun.b1;     % position of nucleaus in z
data{5} = cfun.b1-cfun.c1;  % lower end of nucleaus (z coord.; mean-SD)
data{6} = cfun.b1+cfun.c1;     % upper end of nucleau (mean+SD)
data{7} = nucleus_length; % length of nucleus (4 sigma)
data{8} = nucleus_depth; % depth of nucleus

%keyboard

end