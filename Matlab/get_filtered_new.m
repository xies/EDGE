function [cellsf,filt] = get_filtered_new(cells,ext,llp,hhp,beta)

% Bandpass-filtering of image. In Fourier Domain. Returns filtered image
% and filter kernel. 

a=size(cells);

bg=zeros(ext);     % put on square array of size 2^N
bg(floor((ext-a(1))/2) + 1:floor((ext-a(1))/2)+a(1),floor((ext-a(2))/2)+1:floor((ext-a(2))/2)+a(2))=cells;

[xx,yy] = meshgrid(1:ext,1:ext);   % koordinate system
xx=xx-ext/2;
yy=yy-ext/2;

% Fermi bandpass filter
lp=ext/llp;
filt=1./(1+exp((sqrt(xx.^2+yy.^2)-lp)./(beta)));

hp=ext/hhp;
filt=filt-(hhp~=0)*1./(1+exp((sqrt(xx.^2+yy.^2)-hp)./(beta)));

% filtering in Fourier domain
cellsf=real(ifft2(fft2(bg).*fftshift(filt)));
cellsf=cellsf((ext-a(1))/2+1:(ext-a(1))/2+a(1),(ext-a(2))/2+1:(ext-a(2))/2+a(2));
cellsf=cellsf-(hhp~=0)*mean(mean(cellsf));
filt=filt((ext-a(1))/2+1:(ext-a(1))/2+a(1),(ext-a(2))/2+1:(ext-a(2))/2+a(2));
