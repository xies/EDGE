clear all; close all; clc;

%%
% data_set = 'Histone032';
data_set = '090309 MembCherry HistoneGFP';
cell_inds = []; % all cells
measurement = 'Membranes--vertices--# of neighbors';
layers_from_top = 3; 
masterlayer = 9;

%%
% '1-15-2011, 3;29 PM'

[data_all data_apical data_basal data_middle] = ...
    extract_measurement(data_set, measurement, cell_inds, layers_from_top, masterlayer);

%% take distribution over cells
% out = [];
% %% 
% out(:,1)=neighbor_results(data_apical, 'apical');
% %%
% out(:,2)=neighbor_results(data_basal,  'basal');
% %%
% out(:,3)=neighbor_results(data_middle, 'middle');


%% matthias wants the image thing. OK
plot_neighbor_image(diff(data_apical), 5);
title('nneighbor diff vs. time, apical');
%%
d = diff(data_middle);
d = d(5:end-10,:);  % cut off ends, very liberally...
% for i = 1:size(d,2)  % cells
%     last = 0;
%     for j = 1:size(d,1)  % times
%         if abs(d(j,i)) == 1  % if change
%            if last
%                d(j,i)=0;
%            elseif abs(d(min(j+1,size(d,1)),i)) == 1  % the next one
%                d(j,i)=0;
%            end
%            last = 1;
%         else
%             last = 0;
%         end
%     end
% end

plot_neighbor_image(d, 17);
title('nneighbor diff vs. time, middle');
%%
plot_neighbor_image(diff(data_basal),  7);
title('nneighbor diff vs. time, basal');

%% 
plot_neighbor_image(diff(squeeze(data_all(1, :, :))), 102)
title('nneighbor diff vs z, time 10');
%%

% for cells that have a transition in z, plot the transition location vs
% time

data2 = data_all(:, 1:end-3,:); % get rid of junk at the top

% my_sum(diff(squeeze(data_all(1, :, :)))')
%%
c=138;  %158 102 141 161 (148, 151, 127 at time 10)  % 158, 151,127 nicest  131 ok
% 66, 131 at time 20 138 143 150 (106)
% time 0: 155, 112, 66, 138
d = diff(data2(:,:,c)');
[z1 t1] = find(d < 0);
[z2 t2] = find(d > 0);
% z(6)=[]; t(6) =[];
figure; plot(t1, z1, '.b'); hold on; plot(t2, z2, '.r');
legend('apical->basal vertex gains', 'apical-> basal vertex losses')
xlabel('time')
ylabel('depth of vertex change')
title(['depth of neighbor number change vs time for cell' num2str(c)]);
% vertex-exchange seems to move downwards... and all seem to be RED
% red means d > 0, means you GAIN a vertex moving up, or lose one as you
% move basally