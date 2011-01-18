clear all;  clc; %close all;

%% look at the spatial patterning of nucleus drops


data_set = '090309 MembCherry HistoneGFP';
cell_inds = '1-17-2011, 1;22 PM'; % all cells
cell_inds = [];
measurement = 'Nuclei--nuclei_intensity--nuclei intensity';
layers_from_top = 3;


[data_all data_apical data_basal data_middle neighborhood cell_inds] = ...
    extract_measurement(data_set, measurement, cell_inds, layers_from_top);

%%


tlength = size(data_all, 1);
zlength = size(data_all, 2);
clength = size(data_all, 3);

nuclear_position = NaN(tlength,clength);
nucleus_drop = NaN(clength, 1);
for j = 1:clength
    for i = 1:tlength
        nuclear_position(i,j) = get_nuclear_position(data_all, i, j);
    end
    nucleus_drop(j) = find_nucleus_drop(nuclear_position(:,j));
end

%% do the correlation thing with dropping times
norder = 5;
% neigh = extract_neighborhood(data_set, cell_inds, norder);
load neighbor_info_090309

%%
masterlayer = 9;
mastertime = 1;
dropping_dt = NaN(length(cell_inds), norder);
dropping_xcor = NaN(length(cell_inds), norder, 2*tlength-1);

nuclear_position_noNaN = nuclear_position;
nuclear_position_noNaN(isnan(nuclear_position_noNaN)) = 0;

for i = 1:length(cell_inds)
    droptime = nucleus_drop(i);
    
    % compute dropping_dt
    if ~isnan(droptime)
        for j = 1:norder
            dropping_dt(i, j) = ...
                my_mean( abs( droptime - nucleus_drop( neigh{i, droptime, masterlayer, j} ) ).' );
        end
    end   
    
    % compute xcorr
    for j = 1:norder
        cor = zeros(2*tlength-1, 1);
        for k = 1:length(neigh{i, mastertime, masterlayer, j})
            
            newcorr = xcorr(nuclear_position_noNaN(:,i), ...
                      nuclear_position_noNaN(:,(neigh{i, mastertime, masterlayer, j}(k))), 'coeff');
            cor = cor + newcorr;
        end
        dropping_xcor(i, j, :) = cor / length(neigh{i, mastertime, masterlayer, j});
    end
    
end

dropping_dt_avg = my_mean(dropping_dt.');

%% cross corr of nucleus position
dropping_xcor_avg = NaN(norder, size(dropping_xcor, 3));
for i = 1:size(dropping_xcor, 3)
    dropping_xcor_avg(:, i) = my_mean(dropping_xcor(:,:,i)');
end
figure; plot(-tlength+1:tlength-1, dropping_xcor_avg');
legend('1','2','3','4','5');
xlabel('delay')
ylabel('normalized cross correlation')
title('cross correlation of nucleus position over time for different neighbor orders')

%% cross corr at zero lag vs order
figure; plot(1:norder, dropping_xcor_avg(:, tlength-1), '.');
xlabel('neighbor order');
ylabel('cross correlation at zero lag')

%% avg bar graph of delta t vs neighbor order
figure; bar(1:norder, dropping_dt_avg*19.5/60)
xlabel('neighbor order')
ylabel('average absolute value of time difference (minutes');

%% movie of bar graphs with delta t vs neighbor order
for i = 1:size(dropping, 1);
    if any(isnan(dropping(i,:)))
        continue;
    end
    figure(111)
    bar(1:norder, dropping(i,:)');
    title(i);
    pause(0.3)
end

%% x-corr with velocity
%%
masterlayer = 9;
mastertime = 1;
dropping_dt = NaN(length(cell_inds), norder);
dropping_v_xcor = NaN(length(cell_inds), norder, 2*tlength-1-2);


nuclear_v = diff(nuclear_position);
nuclear_v_noNaN(isnan(nuclear_v)) = 0;


for i = 1:length(cell_inds)
    % compute xcorr
    for j = 1:norder
        cor = zeros(2*tlength-1-2, 1);
        for k = 1:length(neigh{i, mastertime, masterlayer, j})
            
            newcorr = xcorr(nuclear_v_noNaN(:,i), ...
                      nuclear_v_noNaN(:,(neigh{i, mastertime, masterlayer, j}(k))), 'coeff');
            cor = cor + newcorr;
        end
        dropping_v_xcor(i, j, :) = cor / length(neigh{i, mastertime, masterlayer, j});
    end
    
end

%% cross corr of neucleus velocity
dropping_v_xcor_avg = NaN(norder, size(dropping_v_xcor, 3));
for i = 1:size(dropping_v_xcor, 3)
    dropping_v_xcor_avg(:, i) = my_mean(dropping_v_xcor(:,:,i)');
end
figure; plot(-tlength+1+1:tlength-1-1, dropping_v_xcor_avg');
legend('1','2','3','4','5');
xlabel('delay')
ylabel('normalized cross correlation')
title('cross correlation of nucleus velocity over time for different neighbor orders')





%% graphics junk below

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
    pause(1);
    
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
        tosort = [x(maxi) maxi(:)];
        tosort = sortrows(tosort);
        maxi = tosort(end,2);
    end
    
    if ~isempty(maxi)
        % just a sanity check: the peak should be near the top!
        if x(maxi) < my_mean(x') % or median probably better
            maxi = [];
        end
        if sum(~isnan(x)) < 5  % want at least 5 good points to get an estimate
            maxi = [];
        end
    end
    
    plot(maxi, x(maxi), '.r');
    hold off
    title(num2str(i))
    pause(0.5);
end

if isempty(maxi)
    position = NaN;
else
    position = maxi;
end