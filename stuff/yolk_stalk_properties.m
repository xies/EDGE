function [data names units] = yolk_stalk_properties(embryo, getYolk, t, z, c, dx, dz, dt)
% extracts information about yolk stalks at (t, z, c)
% uses functions ring.m, projimg.m
 

% the names
names{1} = 'intensity';
names{2} = 'x-positions';
names{3} = 'y-positions';
names{4} = 'anisotropy';
names{5} = 'angle';
names{6} = 'volume outside cell';
names{7} = 'hole diameter';
names{8} = 'largest hole diameter';

% the units
units{1} = 'arbitrary';
units{2} = 'microns';
units{3} = 'microns';
units{4} = '';
units{5} = 'rad';
units{6} = 'microns^3';
units{7} = 'microns';
units{8} = 'microns';


% for myosin projection, only the top layer needs to be processed
if c ~= 410
% if (c <= 358 || c>369)
     data = num2cell(NaN(8, 1));
     return;
 end



%% do for every layer
yolk = getYolk(t, z);
[cell_img R] = drawCellSmall(embryo, t, z, c); % draw the Cell

yolk = yolk(R(1):R(1)+size(cell_img, 1)-1, R(2):R(2)+size(cell_img, 2)-1);
yolk_intensity = sum(sum(yolk(cell_img)));

data{1} = yolk_intensity;
data{2} = NaN;
data{3} = NaN;
data{4} = NaN;
data{5} = NaN;
data{6} = NaN;
data{7} = NaN;


%% whole-cell properties

% get layers that were tracked; all others are NaN in x_values 
centroids = Cell.centroidStack(embryo.getCellStack(c, t)) * dx;
x_values = centroids(:, 2);
tracked = ~isnan(x_values);

if (z == embryo.masterLayer && sum(~isnan(x_values))>20)

    zlt = embryo.lowestTracked(c, t); % lowest tracked
    ind_tracked = find(tracked);
    zslt = ind_tracked(2)-1;   % second lowest tracked
    nup = 4;                    % select cell outline nup above lowest tracked 
    znlt = ind_tracked(1+nup)-1;   % n lowest tracked
    
    % for yolk stalk, only the layers below "nabv" above and 
    % "nbel" below the lowest layer will be processed        
    nabv=5;
    nbel=7;
    if (embryo.translateZ(zlt)+1<nbel+1)
        nbel=embryo.translateZ(zlt);
    end
    co_out = zeros(8,nabv+nbel+1);
    
    %% get tilt at lower end of cells
    n_sel = round(6/dz);   % number of z layers chosen to calculate tilt
    
    % extract coordinates of centroids
    centroids = Cell.centroidStack(embryo.getCellStack(c, t)) * dx;
    y_values = centroids(:, 1);
    x_values = centroids(:, 2);
    z_values = (1:length(x_values)) * dz;
    z_values = z_values(:);
    
    % select coordinates of 'n_sel' closest layers to z that are not NaN
    ind_an = find(~isnan(x_values));
    z_dist = abs(z_values(embryo.translateZ(zlt)+1)-z_values(ind_an));
    [q,ind] = sort(z_dist);
    ind = ind(1:min(n_sel,length(ind)));   % use 'n_sel' layers or less
    x_sel = x_values(ind_an(ind));
    y_sel = y_values(ind_an(ind));
    z_sel = z_values(ind_an(ind));
    
    % compute tilt via principle components
    X = [x_sel,y_sel,z_sel];
    [coeff,score,roots] = princomp(X);
    
    % first component, vector pointing towards largest variation of centroids
    v = coeff(:,1);
    if (v(3)<0)   % and towards the top of the cell
        v=-v;
    end
    
    % fitted x and y centroid coordinates (in pixel) based on tilt
    meanX = mean(X,1);
    zbase = z_values(embryo.translateZ(zlt)+1-nbel:...
                     embryo.translateZ(zlt)+1+nabv);
    z_offsets = meanX(1,3)-zbase;     % offset in microns
    xc_fit=meanX(1,1)-z_offsets*v(1)/v(3);
    yc_fit=meanX(1,2)-z_offsets*v(2)/v(3);
    
    ang_tot = atan(sqrt(v(1)^2+v(2)^2)/abs(v(3)));
    ang_dir = atan2(v(2),v(1));

    
 
    %%
    c
    
    for ii=1:nabv+nbel+1
        z_i = zlt-nbel-1+ii;
        
        % get yolk for all cells
        yolk = getYolk(t, z_i);
        if (z_i>znlt && tracked(embryo.translateZ(z_i)+1))
            [roic, R] = drawCellSmall(embryo, t, z_i, c);  % roic is cell_img
            % (R(1) is offset in vertical, R(2) in horizontal direction)
        else
            [roic, R] = drawCellSmall(embryo, t, znlt, c);  % roic is cell_img
            % (R(1) is offset in vertical, R(2) in horizontal direction)
            z_offset = z_values(embryo.translateZ(znlt)+1)-...
                       z_values(embryo.translateZ(z_i)+1);     % offset in microns
            xsh=round(-z_offset*v(1)/v(3)/dx);
            ysh=round(-z_offset*v(2)/v(3)/dx);
            R=R+[ysh xsh];        
        end
        npad = 2;
        roic = padarray(roic,[npad npad]); % increase region around cell
        roic = bwmorph(roic,'thicken',npad);
        R = R-[npad npad];
        yxs = size(roic);
            
        
        % extract yolk inside cell
        yolkc = double(yolk(R(1):R(1)+yxs(1)-1,R(2):R(2)+yxs(2)-1));
        yolkc = yolkc.*roic;
        yolkc_raw = yolkc;
        yolkc = yolkc/mean(mean(yolkc(roic))); % normalize
        
        % project on plane perpendicular to cell axis
        roic = projimg(roic,ang_tot,ang_dir,...
                      [xc_fit(ii)/dx-R(2) yc_fit(ii)/dx-R(1)]);
        yolkc = projimg(yolkc,ang_tot,ang_dir,...
                      [xc_fit(ii)/dx-R(2) yc_fit(ii)/dx-R(1)]);
        yolkc = yolkc.*roic/cos(ang_tot);
keyboard
        
        %initial paramters for fitting:
        ext=(2^ceil(max(log2(size(yolkc)))));       %smooth yolkin inside cell
        [yolkf,filt] = get_filtered_gauss(yolkc,ext,2,0);
        [roicf,filt] = get_filtered_gauss(roic,ext,2,0);
        yolkf(roic>0)=yolkf(roic>0)./roicf(roic>0);
        yolkf=yolkf.*double(roic);
        
        yolk_th=yolkf>0.8;                    %and threshold
        yolki = bwlabel(yolk_th,4);
        yolki(yolki>1)=1;
        
        propty = regionprops(yolki,'Centroid');     %get initial parameters for fitting
        cents = [propty.Centroid];
        nspots = size(cents,2)/2;
        cents = transpose(reshape(cents,[2 nspots]));
        
        
        propty = regionprops(yolki,'MajorAxisLength');
        mxax = transpose([propty.MajorAxisLength]);
        propty = regionprops(yolki,'MinorAxisLength');
        mnax = transpose([propty.MinorAxisLength]);
        propty = regionprops(yolki,'Orientation');
        orient = transpose([propty.Orientation]);
        orient = -orient/90*pi/2;
        
        
        %     ng=length(mxax);
        
        xe0=1;      % set paramters and upper and lower bounds
        ub=3;
        lb=0;
        
        %xe0=[xe0 1 cents(1,1) cents(1,2) ...
        cenx = xc_fit(ii)/dx-R(2); 
        ceny = yc_fit(ii)/dx-R(1);
        xe0=[xe0 1 cenx ceny...
            max([mxax(1)/3 1/dx]) ...
            max([mnax(1)/3 1/dx])/max([mxax(1)/3 1/dx])...
            orient 2/dx 0.1 0.2];
        maxradius=max([yxs(2) yxs(1)]);
        ub=[ub 6 cenx+1/dx ceny+1/dx maxradius/2 1.4 300 1.5/dx 1 300];
        lb=[lb 0 cenx-1/dx ceny-1/dx 0 0.7 -300 0 0 -300];
        
        
        if ~isempty(orient)
            opts = optimset('display', 'off');
            x = lsqcurvefit(@ring, xe0, roic, yolkc, lb, ub, opts);   %fit
            cf=x(2:end);
            cf(5)=cf(4)*cf(5);
            co_out(1,ii) = dx*(R(2)+cf(2)); % 'Stalk x position'
            co_out(2,ii) = dx*(R(1)+cf(3)); % 'Stalk y position'
            if (cf(4)>=cf(5))
                co_out(3,ii) = dx*cf(4);      % 'Stalk larger half axis'
                co_out(4,ii) = dx*cf(5);      % 'Stalk smaller half axis'
                qangle=exp(complex(0,1)*(cf(6)+pi/2));   % shift to angle of smaller axis 
                qangle=atan(imag(qangle)/real(qangle)); % bring angle between -pi/2 and pi/2
            else
                co_out(3,ii) = dx*cf(5);      % 'Stalk larger half axis'
                co_out(4,ii) = dx*cf(4);      % 'Stalk smaller half axis'
                qangle=exp(complex(0,1)*(cf(6)));   % use angle ofsmaller axis
                qangle=atan(imag(qangle)/real(qangle)); % bring angle between -pi/2 and pi/2
            end
            co_out(5,ii) = qangle;         % 'Stalk angle'
            co_out(6,ii) = dx*cf(7);          % 'Stalk width'
            co_out(7,ii) = 2*dx*(cf(5)-cf(7)); % 'Hole diameter'
            yolkc_fit = ring(x,roic);
            %intensity = sum(sum(yolkc_fit>1))/sum(sum(roic));
            intensity = mean(mean(yolkc_raw(roic)));
            co_out(8,ii) = intensity;          % 'Stalk intensity'
        else
            co_out(:,ii)=NaN(7,1);
        end        
    %keyboard
    end
    
    %% calculate properties that characterize yolk stalk
    
    min_ysint_fac = 0.8; %intensity factor below which images are discarded 
    
    % x and y position of center (of highest intensity) of yolk stalk
    [mxint,ind_mx_ys] = max(co_out(8,:));
    data{2} = co_out(1,ind_mx_ys);
    data{3} = co_out(2,ind_mx_ys);

    % angle
    data{4} = co_out(3,ind_mx_ys)/co_out(4,ind_mx_ys);
    data{5} = co_out(5,ind_mx_ys);

    % yolk stalk volume outside of cell
    ys_vol = 0;
    if (nbel>0)
        for ii=1:nbel
            if (co_out(8,ii)>min_ysint_fac*co_out(8,ind_mx_ys))
                ys_ar = pi*(co_out(3,ii)+co_out(6,ii))*...
                           (co_out(4,ii)+co_out(6,ii));
                ys_vol = ys_vol + dz*ys_ar;
            end               
        end    
        data{6}=ys_vol;
    else
        data{6}=NaN;
    end
    
    % smallest hole in yolk stalk
    sm_radius=zeros(nbel+1+nabv,1);
    for ii=1:nbel+1+nabv
        if (co_out(8,ii)>min_ysint_fac*co_out(8,ind_mx_ys))
            sm_radius(ii) = co_out(4,ii);
        end
    end
    data{7}=min(sm_radius(sm_radius~=0));    
    data{8}=max(sm_radius(sm_radius~=0));    
    
    
    %keyboard    
    

end

end

%%

function F = ring(x,xdata)
% generates function used to fit the yolk stalk: 
% constant plus an anisotropic rotated ring


a=size(xdata);
[xk,yk] = meshgrid(1:a(2),1:a(1));


F=xk*0+x(1);

height=x(2);        % height
cx=x(3);            % x-position of ring
cy=x(4);            % y-position of ring
radius=x(5);        % radius of ring in x direction
facy=x(6);          % factor in front of radius of ring in y direction 
theta=x(7);         % rotation angel of ring
sig=x(8);           % width of ring
dampf=x(9);         % factor for damping one side of ring (ramp)
angdamp=x(10);      % in direction of angle


xkt=xk-cx;
ykt=yk-cy;
xkturn=xkt*cos(theta)+ykt*sin(theta);
ykturn=-xkt*sin(theta)+ykt*cos(theta);
ykturn=ykturn/facy;

% create ring
g=height*exp(-(sqrt(xkturn.^2+ykturn.^2)-radius).^2/2/sig^2);

% ramp damping
damp=xkt*cos(angdamp)+ykt*sin(angdamp);
damp=damp/max(abs(damp(:)));
damp=damp*0+1+damp*dampf;

% add background ring and ramp damping
F=F+g.*damp;

F=F+(~xdata(round(cy),round(cx)))*10;    % penality for centering Gaussians on background
F=F.*double(xdata);

end


%%
function mcimg = projimg(cimg,ang_tot,ang_dir,cents)
% projets cell image (cimg) on plane perpendicular to cell axis
% ang_tot is angle of cell axis with z axis
% ang_dir is angle of cell axis in xyplane 
% 

a = size(cimg);   % size of image

omega = [cos(ang_dir)  sin(ang_dir) ; ...  % rotates to
        -sin(ang_dir)  cos(ang_dir)];       % x axis

W = [cos(ang_tot) 0; 0 1];   % stretch matrix (along x axis)     
     
W=omega'*W*omega;            % and back
W = [W; 0 0];                 % adds movement

    
T = maketform('affine',W);   % build the transform and perform

[mcimg, xdata, ydata] = imtransform(cimg,T,...
         'udata',[-cents(1)+1 -cents(1)+a(2)],...
         'vdata',[-cents(2)+1 -cents(2)+a(1)],...
         'xdata',[-cents(1)+1 -cents(1)+a(2)],...
         'ydata',[-cents(2)+1 -cents(2)+a(1)]);



end


