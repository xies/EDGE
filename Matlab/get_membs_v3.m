function cellsi = get_membs_v3(cells, lp, hp, th)

% [Ys Xs] = size(cells);

% cells = scale(cells);

% preprocessing
cells = cells - mean(mean(cells));     %subtract mean  

%truncate values at +- trunc_th standard deviations
sigma = std(cells(:));
trunc_th = 2;

% figure, imshow(cells);

cells(cells>=trunc_th*sigma)=trunc_th*sigma;  %find wherever cells is greater than this threshold and set it to this threshold
cells(cells<=-trunc_th*sigma)=-trunc_th*sigma;  %same for below the bottom threshold
%filtering


% figure, imshow(cells);

[cellsf,filt] = get_filtered(cells,lp,hp);  %cellsf is the filtered cells array
%    figure, imshow(scale(cellsf))


% thresholding and binary and start skeletonizing
bw_thresh = th * std(cellsf(:));

cellsb = im2bw(cellsf - bw_thresh, 0); 


%cellsb = bwmorph(cellsb,'close');

    %figure, imshow(scale(cellsb))

  %  keyboard;
%cellsX = watershed(cellsb);
cellsX = cellsb;
%cellsX = - cellsX+1;



    %figure, imshow(scale(cellsb))
%keyboard;

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
ddx=5;ddy=5; 
mask(1:ddy,:) = 0; mask(Ys-ddy:Ys,:) = 0; mask(:,1:ddx) = 0; mask(:,Xs-ddx:Xs) = 0;

cellsX(mask==0)=1;

cellsi = cellsX == 1;

% skeletonize
cellsi = bwmorph(cellsi,'shrink',Inf);  %OFF    %thins to 1 px
%    figure, imshow(cellsi);

cellsi = bwmorph(cellsi,'clean');   %gets rid of single dots



% label cells by integers; every cell gets addressed
cellsi = -cellsi+1;  %flipping ones and zeros
cellsi = logical(cellsi);

cellsi(mask==0)=1;

cellsi = ~cellsi;
