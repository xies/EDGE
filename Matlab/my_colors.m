function col = my_colors(x)
% here i define what colors i want to use for the multiple channels and
% also for selecting a large number of cells. i want to avoid blue because
% i already use this for the polygons. all other colors are OK (i don't
% care much about the borders)

% this function should be defined for any positive integer, and is
% deterministic (important!)

% blue is [0 0 1] --> avoid this
% other than this criterion, i pretty much just picked randomly

% this function can handle an array if x-values and outputs an n x 3 array
% where n = length(x)


% R = sin(x) + .8*cos(x);
% G = (pi*sqrt(x) - x)/10 + 0.05*cos(x/7.7);
% B = log2(x+1)/24 + sin(x/5) - 0.3;
% 
% %%% set the first few by default
% % if x = 1, pick red
% R(x == 1) = 1;
% G(x == 1) = 0;
% B(x == 1) = 0;
% % if x = 2, pick green
% R(x == 2) = 0;
% G(x == 2) = 1;
% B(x == 2) = 0;
% % if x = 3, pick some weird purpose
% R(x == 3) = 0.45;
% G(x == 3) = 0.15;
% B(x == 3) = 0.4;
% 
% 
% 
% col = [R(:) G(:) B(:)];
% 
% for i = 1:size(col, 1)
%     col(i, :) = col(i, :) - min(0, min(col(i, :)));  % make sure they are non-negative
%     col(i, :) = col(i, :) / max(col(i, :));  % make sure they are all less than 1
% end

% H = (x.^(-0.2))*cos((x-1)/1.3).^2 + cos((x-1)/19.5)/20;
H = abs(3-mod(x-1, 6))/3 + sin((x-1)/2.23)^2/7;
S = ones(length(x), 1);
V = ones(length(x), 1);

H = min(1, max(0, H));  % force it to be on [0 1]

[R G B] = hsv2rgb([H(:) S(:) V(:)]);

col = [R B G];  % avoid blue at first