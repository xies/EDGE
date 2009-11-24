function [cellsf,filt] = get_filtered(cells,llp,hhp)

%cells = array of cells
%ext = size of background array
%llp = low pass cutoff (pixels)
%hhp = high pass cutoff (pixels)
%beta = sharpness of filter

% Bandpass-filtering of image. In Fourier Domain. Returns filtered image
% and filter kernel. 

a = size(cells);


beta = 8;   % sharpness of filter
ext = 2^(ceil(log2(max(a)))); % smallest power of 2 bigger than Xs and Ys


bg=zeros(ext);     % put on square array of size 2^N

bg(floor((ext-a(1))/2) + 1:floor((ext-a(1))/2)+a(1), floor((ext-a(2))/2)+1:floor((ext-a(2))/2)+a(2))=cells;

[xx,yy] = meshgrid(1:ext,1:ext);   % coordinate system
xx=xx-ext/2;
yy=yy-ext/2;

% Fermi bandpass filter
lp=ext/llp;
hp=ext/hhp;
filt=1./(1+exp((sqrt(xx.^2+yy.^2)-lp)./(beta)));
filt=filt-1./(1+exp((sqrt(xx.^2+yy.^2)-hp)./(beta)));


%figure, imshow(fftshift(filt));

% multiplying by filter in Fourier domain
fourier_temp = fft2(bg);

%figure, imagesc(abs(fourier_temp) - mean(mean(abs(fourier_temp))));

fourier_temp = fourier_temp.*fftshift(filt);

%figure, imshow(abs(fourier_temp));

cellsf=real(ifft2(fourier_temp));   %filt is the circle thing in figure 1

cellsf=cellsf(floor((ext-a(1))/2)+1:floor((ext-a(1))/2)+a(1), ...
    floor((ext-a(2))/2)+1:floor((ext-a(2))/2)+a(2));
cellsf=cellsf-mean(mean(cellsf));
filt=filt(floor((ext-a(1))/2)+1:floor((ext-a(1))/2)+a(1),floor((ext-a(2))/2)+1:floor((ext-a(2))/2)+a(2));

% % pad with zeros if necessary to get an even number
% if mod(size(cells, 1), 2)==1
%     cellsf = [cellsf; zeros(1, size(cells, 2))];
% end
% if mod(size(cells, 2), 2)==1
%     cellsf = [cellsf zeros(size(cells, 1), 1)];
% end