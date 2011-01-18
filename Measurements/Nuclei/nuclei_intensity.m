function [data names units] = nuclei_intensity(embryo, getNuclei, t, z, c, dx, dz, dt, other)

% the names
names{1} = 'nuclei intensity';
names{2} = 'nuclear position';

% the units
units{1} = 'intensity units';
units{2} = 'microns';


% nuclei = getNuclei(t, z);
% nuclear_intensity = sum(sum(nuclei(cell_img))); 

if z == embryo.masterLayer
    nuclear_intensity_profile = zeros(embryo.z, 1);
    for z_i = embryo.lowestTracked(c, t) : ...
            my_sign(embryo.highestTracked(c, t) - embryo.lowestTracked(c, t)) : ...
            embryo.highestTracked(c, t)
        
        % draw the Cell
        [cell_img R] = drawCellSmall(embryo, t, z, c);
        
        nuclei = getNuclei(t, z_i);
        nuclei = nuclei(R(1):R(1)+size(cell_img, 1)-1, R(2):R(2)+size(cell_img, 2)-1);
        nuclear_intensity = sum(sum(nuclei(cell_img))); 
        nuclear_intensity_profile(embryo.translateZ(z_i)+1) = nuclear_intensity; 
        
    end

    [~, max_ind] = max(nuclear_intensity_profile);
    % max_ind = embryo.unTranslateZ(max_ind); 

    data{1} = nuclear_intensity_profile(embryo.translateZ(z)+1);
    data{2} = max_ind * dz;
else
    % draw the Cell
    [cell_img R] = drawCellSmall(embryo, t, z, c);
    
    nuclei = getNuclei(t, z);
    nuclei = nuclei(R(1):R(1)+size(cell_img, 1)-1, R(2):R(2)+size(cell_img, 2)-1);
    nuclear_intensity = sum(sum(nuclei(cell_img))); 
    
    data{1} = nuclear_intensity;
    data{2} = NaN;
end