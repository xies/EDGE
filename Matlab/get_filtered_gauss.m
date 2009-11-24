function [cellsf,kernk] = get_filtered_gauss(cells,ext,slp,shp)

% Bandpass-filtering of image. In Fourier Domain. Returns filtered image
% and filter kernel. 

a=size(cells);

bg=zeros(ext);     % put on square array of size 2^N
bg(floor((ext-a(1))/2) + 1:floor((ext-a(1))/2)+a(1),floor((ext-a(2))/2)+1:floor((ext-a(2))/2)+a(2))=cells;

[xx,yy] = meshgrid(1:ext,1:ext);   % koordinate system
xx=xx-ext/2-1;
yy=yy-ext/2-1;

% Gauss bandpass filter
kernk=xx*0+1;
if (shp~=0)
    kern=1/(2*pi*shp^2)*exp(-(xx.^2+yy.^2)/2/shp^2);
    kernk=kernk-real(fft2(fftshift(kern)));
end
if (slp~=0)
    kern=1/(2*pi*slp^2)*exp(-(xx.^2+yy.^2)/2/slp^2);
    kernk=kernk.*real(fft2(fftshift(kern)));
end 
    
% filtering in Fourier domain
cellsf=real(ifft2(fft2(bg).*kernk));
cellsf=cellsf(floor((ext-a(1))/2) + 1:floor((ext-a(1))/2)+a(1),floor((ext-a(2))/2)+1:floor((ext-a(2))/2)+a(2));
kernk=fftshift(kernk);
kernk=kernk(floor((ext-a(1))/2) + 1:floor((ext-a(1))/2)+a(1),floor((ext-a(2))/2)+1:floor((ext-a(2))/2)+a(2));
