function cellsi = get_membs_v3(cells, lp, hp, th)

% preprocessing
cells = cells - mean(mean(cells));     %subtract mean  

%truncate values at +- trunc_th standard deviations
sigma = std(cells(:));
trunc_th = 2;

cells(cells>=trunc_th*sigma)=trunc_th*sigma;  %find wherever cells is greater than this threshold and set it to this threshold
cells(cells<=-trunc_th*sigma)=-trunc_th*sigma;  %same for below the bottom threshold
%filtering


[cellsf,filt] = get_filtered(cells,lp,hp);  %cellsf is the filtered cells array

% thresholding and binary and start skeletonizing
bw_thresh = th * std(cellsf(:));

cellsb = im2bw(cellsf - bw_thresh, 0); 

%cellsX = watershed(cellsb);
cellsX = cellsb;




%cellsX = - cellsX+1;

% first, confine to ROI, i.e.region within boundary, by multpl. mask
% mask = bwmorph(cellsb,'thin',Inf);     % prepare mask
% mask = bwmorph(mask,'shrink',Inf);     % ...
% mask = bwmorph(mask,'clean');

% figure, imshow(mask);
% mask = -mask + 1;
% mask = logical(mask);
% mask = bwlabel(mask,4);
% mask((mask==1)) = -1;
% mask((mask>=0)) = 1;
% mask((mask==-1)) = 0;        % mask 1 inside ROI, else 0

%mask=ones(Ys,Xs);
%ddx=5;ddy=5; 
%mask(1:ddy,:) = 0; mask(Ys-ddy:Ys,:) = 0; mask(:,1:ddx) = 0; mask(:,Xs-ddx:Xs) = 0;
% shrink ROI by multiplying mask with shifted version

%mask=mask.*circshift(mask,[ddy,ddx]).*circshift(mask,[-ddy,-ddx]);
%   %ORIGINAL MASK (above)

%cellsX=cellsX.*mask;

[Ys Xs] = size(cellsX);  % in case get_filtered returns something
                         % different in size by 1 pixel

mask=ones(Ys,Xs);
ddx=2;ddy=2; 
mask(1:ddy,:) = 0; mask(Ys-ddy+1:Ys,:) = 0; mask(:,1:ddx) = 0; mask(:,Xs-ddx+1:Xs) = 0;
%keyboard

cellsX(mask==0)=1;

cellsi = cellsX == 1;

% skeletonize
cellsi = bwmorph(cellsi,'shrink',Inf);    %thins to 1 px

cellsi = bwmorph(cellsi,'clean');   %gets rid of single dots


cellsi(mask==0)=0;
% label cells by integers; every cell gets addressed
cellsi = logical(cellsi);
