% Get the membranes from a raw image of cells. First filter (lp, hp), then
% threshold (th), then do all the skeletonization and cleaning up. 

function cellsi = get_membs_v3(cells, lp, hp, th)


% preprocessing
cells = cells - mean(mean(cells));     %subtract mean  

%truncate values at +- trunc_th standard deviations
sigma = std(cells(:));
trunc_th = 2;

cells(cells>=trunc_th*sigma)=trunc_th*sigma;  %find wherever cells is greater than this threshold and set it to this threshold
cells(cells<=-trunc_th*sigma)=-trunc_th*sigma;  %same for below the bottom threshold
%filtering


% relatively slow?
[cellsf,filt] = get_filtered(cells,lp,hp);  %cellsf is the filtered cells array



% thresholding and binary and start skeletonizing
bw_thresh = th * std(cellsf(:));

cellsX = im2bw(cellsf - bw_thresh, 0); 

[Ys Xs] = size(cellsX);  % in case get_filtered returns something
                         % different in size by 1 pixel
% throw away the outer edges
mask=ones(Ys,Xs);
ddx=2;ddy=2; 
mask(1:ddy,:) = 0; mask(Ys-ddy+1:Ys,:) = 0; mask(:,1:ddx) = 0; mask(:,Xs-ddx+1:Xs) = 0;

cellsX(mask==0)=1;
cellsi = cellsX == 1;

% skeletonize (relatively slow?)
cellsi = bwmorph(cellsi,'shrink',Inf);    %thins to 1 px
cellsi = bwmorph(cellsi,'clean');   %gets rid of single dots

cellsi(mask==0)=0;
% label cells by integers; every cell gets addressed
cellsi = logical(cellsi);