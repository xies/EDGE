function [smo_out,hf_out,roc_out]=smooth_rochange(data,sigma,dt)
% input
% data= 2D input array (n timeseries), smoothed along FIRST dimension
% sigma= width of Gaussian smoothing kernel
% dt= time resolution, time per pixel
% output
% smo= smoothed vertion of q
% hf= high freq part of q (q-smo)
% roc= rate of change (calc. based on smo)

smo_out = zeros(size(data));
hf_out = zeros(size(data));
roc_out = zeros(size(data));

for I = 1:size(data, 2)
    q = data(:, I).';
    q_old = q;  
    firstNaN = find(isnan(q), 1);
    if (isempty(firstNaN))
        firstNaN = length(q) + 1;
    end
    % after the first NaN entry, all data is ignored
    % this is not so ideal
    q(firstNaN:end) = [];
    q_old(firstNaN:end) = NaN;
    %q_old is used to put the smoothed stuff
    % back in the right sized array and then pacakge it
    if firstNaN == 1
        continue;
    end
  %%  
    % create  backgraound array
    framesize=fix(4*sigma);
    aa = size(q)+[0 2*framesize];
    [xx,yy] = meshgrid(1:aa(2),1:aa(1)); 
    xx=xx-aa(2)/2;

    % filter kernel
    filt=exp(-xx.^2/2.0/sigma^2)/sqrt(2*pi)/sigma;
    filt=circshift(filt,[0,fix(aa(2)/2)+1]);
    filtk=fft(filt,[],2);

    if (sigma==0)
       filtk=ones(aa(1),aa(2)); 
    end

    %norm to avoid edge effects
    norm=zeros(aa(1),aa(2));
    norm(:,framesize+1:end-framesize)=1;
    norm=real(ifft(filtk.*fft(norm,[],2),[],2));
    norm=norm(:,framesize+1:end-framesize);

    %smoothing
    bg=zeros(aa(1),aa(2));
    bg(:,framesize+1:end-framesize)=q;
    fbg = fft(bg,[],2);
    bgsm=real(ifft(fbg.*filtk,[],2));
    smo=bgsm(:,framesize+1:end-framesize)./norm;
    hf = q-smo;
    if (sigma==0)
       hf=hf*0; 
    end
%%

    % rate-of-change
    ll=circshift(smo,[0 -1]);
    rr=circshift(smo,[0 1]);
    roc=(ll-rr)/dt/2;
    roc(:,1)=NaN;
    roc(:,end)=NaN;
    
    z = q_old;
    z(1:firstNaN-1) = smo;
    smo_out(:, I) = z.';
    z = q_old;
    z(1:firstNaN-1) = hf;
    hf_out(:, I)  = z.';
    z = q_old;
    z(1:firstNaN-1) = roc;
    roc_out(:, I) = z.';
    
end

