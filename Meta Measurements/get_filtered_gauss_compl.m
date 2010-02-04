function [cellsf,filt] = get_filtered_gauss(cells,ext,slp,shp)

% Bandpass-filtering of image. In Fourier Domain. Returns filtered image
% and filter kernel. 

a=size(cells);

bg=zeros(ext);     % put on square array of size 2^N
bg((ext-a(1))/2+1:(ext-a(1))/2+a(1),(ext-a(2))/2+1:(ext-a(2))/2+a(2))=cells;

[xx,yy] = meshgrid(1:ext,1:ext);   % koordinate system
xx=xx-ext/2-1;
yy=yy-ext/2-1;

% Gauss bandpass filter
filt=xx*0;
if (slp~=0)
    filt=filt+1/(2*pi*slp^2)*exp(-(xx.^2+yy.^2)/2/slp^2);
else 
    filt(ext/2+1,ext/2+1)=1;
end 
if (shp~=0)
    filt=filt-1/(2*pi*shp^2)*exp(-(xx.^2+yy.^2)/2/shp^2);
end

    
% filtering in Fourier domain
cellsf=(ifft2(fft2(bg).*fft2(fftshift(filt))));
cellsf=cellsf((ext-a(1))/2+1:(ext-a(1))/2+a(1),(ext-a(2))/2+1:(ext-a(2))/2+a(2));
filt=filt((ext-a(1))/2+1:(ext-a(1))/2+a(1),(ext-a(2))/2+1:(ext-a(2))/2+a(2));
