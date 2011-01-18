clear all; close all; clc;

%% look at the spatial patterning of nucleus drops


data_set = '090309 MembCherry HistoneGFP';
cell_inds = '1-17-2011, 1;22 PM'; % all cells
measurement = 'Nuclei--nuclei_intensity--nuclei intensity';
layers_from_top = 3;


[data_all data_apical data_basal data_middle neighborhood] = ...
    extract_measurement(data_set, measurement, cell_inds, layers_from_top);

%%


nuclear_position = zeros(size(data_all,1),size(data_all,3));
nucleus_drop = zeros(size(data_all, 3));
for j = 1:size(data_all,3)
    for i = 1:size(data_all, 1)
        nuclear_position(i,j) = get_nuclear_position(data_all, i, j);
    end
    nucleus_drop(j) = find_nucleus_drop(nuclear_position(:,j));
end

%%




%%

nucleus_drop = zeros(size(data_all, 3));

figure(1);
for i = 1:size(data_all,3)
    nucleus_drop(i) = find_nucleus_drop(nuclear_position(:,i));
    
    nucleus_drop(i)
    i
    
    plot(nuclear_position(:,i));
    hold on;
    if ~isnan(nucleus_drop(i))
        plot(nucleus_drop(i), nuclear_position(nucleus_drop(i),i), '*r');
    end
    hold off;
    title(i);
    pause(1.5);
    
end



%%
c=10; 
figure(1)
for i = 1:size(data_all, 1)
    x = data_all(i, :, c); 
    x=x(:);
    plot(x); 
    hold on
    maxi = find( diff([NaN; x])>0 &  diff([x; NaN])<0   );
    % of these "peaks", find the one at the maximum intensity
    
    if length(maxi) > 1
        tosort = [maxi(:) x(maxi)];
        tosort = sortrows(tosort);
        maxi = tosort(end,1);
    end
    if ~isempty(maxi)
        % just a sanity check: the peak should be near the top!
        if x(maxi) < my_mean(x')
            maxi = [];
        end
        if sum(~isnan(x)) < 5  % want at least 5 good points to get an estimate
            maxi = [];
        end
    end
    
    plot(maxi, x(maxi), '.r');
    hold off
    title(num2str(i))
    pause(.5);
end

if isempty(maxi)
    position = NaN;
else
    position = maxi;
end