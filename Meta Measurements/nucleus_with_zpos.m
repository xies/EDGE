clear all;  clc;
% close all;

%%
data_set = 'Histone032'; masterlayer = 8; dt=5; windowsize = 11;
% data_set = '090309 MembCherry HistoneGFP'; masterlayer=9; dt=19.5; windowsize = 3;
cell_inds = []; % all cells
% measurement = 'Membranes--vertices--# of neighbors';
measurement = 'Nuclei--nuclei_properties--position-z';

layers_from_top = 3;

%% get x y z positions
[data_all,~,~, data_z] = ...
    extract_measurement(data_set, 'Nuclei--nuclei_properties--position-z', cell_inds, layers_from_top, masterlayer);
%%
[~,~,~, data_x] = ...
    extract_measurement(data_set, 'Nuclei--nuclei_properties--position-x', cell_inds, layers_from_top, masterlayer);
[~,~,~, data_y] = ...
    extract_measurement(data_set, 'Nuclei--nuclei_properties--position-y', cell_inds, layers_from_top, masterlayer);
% [~,~,~, data_x] = ...
%     extract_measurement(data_set, 'Membranes--basic_2d--Centroid-x', cell_inds, layers_from_top, masterlayer);
% [~,~,~, data_y] = ...
%     extract_measurement(data_set, 'Membranes--basic_2d--Centroid-y', cell_inds, layers_from_top, masterlayer );


%%

tlength = size(data_z, 1);
clength = size(data_z, 2);
%%


% first, generate the velocities
datavx = NaN(clength, tlength-1);
datavy = NaN(clength, tlength-1);
datavz = NaN(clength, tlength-1);
for i = 1:clength
%     datavx(i,:) = take_derivative(data_x(:,i), windowsize);
%     datavy(i,:) = take_derivative(data_y(:,i), windowsize);
    datavz(i,:) = take_derivative(data_z(:,i), windowsize);
end
%%
% datav = sqrt(datavx.^2 + datavy.^2 + datavz.^2);
datav = datavz;
% datav = sqrt(datavx.^2 + datavy.^2);

%% do the correlation thing
norder = 9;
disp('loading neighbors...')
% neigh = extract_neighborhood_masterlayer(data_set, cell_inds, norder);
% load neighbor_info_090309_ml_9
load neighbor_info_histone032_9
disp('done.');

%% averaging over neighbors FOR EACH cell, THEN averaging over cells
% 
% masterlayer = 9;
% mastertime = 1;
% dropping_xcornolag = NaN(clength, norder+1);
% 
% for i = 1:clength
% i
%     thiscell = data(i,:);
%     thiscell = take_derivative(thiscell);
% %     thiscell = rand(size(thiscell))
%    
%     thiscell = thiscell - my_mean(thiscell);
%     thiscell = thiscell / my_std(thiscell);
%     
%     dot = thiscell.*thiscell;
%     dot = dot(~isnan(dot));
%     corr = sum(dot)/(length(dot)-1);
%     if ~isempty(dot)
%         dropping_xcornolag(i,1) = corr;
%     else
%         dropping_xcornolag(i,:) = NaN(1,norder+1);
%         continue;
%     end
%     
%     % compute xcorr_nolag
%     for j = 1:norder
%         corr = 0;
%         nk = length(neigh{i, mastertime, masterlayer, j});
%         for k = 1:nk % for each cell of that order
%             thatcell = data(neigh{i, mastertime, masterlayer, j}(k),:);
%             
%             thatcell = take_derivative(thatcell);
% %             thatcell = rand(size(thatcell));
%             thatcell = thatcell - my_mean(thatcell);
%             thatcell = thatcell / my_std(thatcell);
% 
%             dot = thiscell.*thatcell;
%             dot = dot(~isnan(dot));
%             if ~isempty(dot)
%                 corr = corr + sum(dot)/(length(dot)-1);
%             end
%         end
%         dropping_xcornolag(i, j+1) = corr/nk;
%     end
%     
% end
% 
% dxc_avg = my_mean(dropping_xcornolag');
% dxc_std = my_std(dropping_xcornolag');

 %% averaging over ALL first order, ALL 2nd order, etc.
% masterlayer = 9;
% mastertime = 1;
% dropping_xcornolag = NaN(norder+1,1);
% 
% dropping_xcornolag(1)=1;
% 
% for j = 1:norder
% 
%     count = 0;
%     corr = 0;
%     for i = 1:clength
%         fprintf('%u %u\n',j, i)
%         thiscell = data(i,:);
%         thiscell = take_derivative(thiscell);
%     %     thiscell = rand(size(thiscell))
% 
%         thiscell = thiscell - my_mean(thiscell);
%         thiscell = thiscell / my_std(thiscell);
% 
% %         dot = thiscell.*thiscell;
% %         dot = dot(~isnan(dot));
% %         corr = sum(dot)/(length(dot)-1);
% %         if ~isempty(dot)
% %             dropping_xcornolag(i,1) = corr;
% %         else
% %             dropping_xcornolag(i,:) = NaN(1,norder+1);
% %             continue;
% %         end
%         
%         nk = length(neigh{i, mastertime, masterlayer, j});
%         for k = 1:nk % for each cell of that order
%             thatcell = data(neigh{i, mastertime, masterlayer, j}(k),:);
%             
%             thatcell = take_derivative(thatcell);
% %             thatcell = rand(size(thatcell));
%             thatcell = thatcell - my_mean(thatcell);
%             thatcell = thatcell / my_std(thatcell);
% 
%             dot = thiscell.*thatcell;
%             dot = dot(~isnan(dot));
%             if ~isempty(dot)
%                 corr = corr + sum(dot)/(length(dot)-1);
%                 count = count + 1;
%             end
%         end
%         
%     end
%     
%     dropping_xcornolag(j+1) = corr/count;
% end
% 
% dxc_avg = dropping_xcornolag;

%%

%% compute xcor separately for each time
% make this big vectors

mastertime = 1;
dropping_xcornolag = NaN(norder+1,tlength-1);

dropping_xcornolag(1,:)=1;
samplesize = zeros(norder+1,1);

for t = 1:tlength-1
    fprintf('time %u of %u\n', t, tlength-1);
    for j = 1:norder

        a=[];b=[];
        for i = 1:clength
            thiscell = datav(i,t);
        %     thiscell = rand(size(thiscell))

            nk = length(neigh{i, t, j});
            for k = 1:nk % for each cell of that order
                thatcell = datav(neigh{i, t, j}(k),t);
                a = [a thiscell];
                b = [b thatcell];
            end

        end
        
        a = a - my_mean(a);
        b = b - my_mean(b);
        a = a / my_std(a);
        b = b / my_std(b);
        dot = a .* b;
        dot = dot(~isnan(dot));
        if isempty(dot)
            dropping_xcornolag(j+1, t) = NaN;
        else
            dropping_xcornolag(j+1, t) = sum(dot)/(length(dot)-1);
        end
    end
end

dxc_avg = dropping_xcornolag;



% figure; plot(dropping_xcornolag.');
% legend('0','1','2','3','4','5');
figure;imagesc((1:tlength-1)*dt/60,0:norder, dropping_xcornolag); colorbar;
ylabel('neighbor order');
xlabel('time (min)')
