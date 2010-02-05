function [data names units] = yolk_stalk_properties(embryo, getYolk, t, z, c, dx, dz, dt, other)
% extracts information about yolk stalks at (t, z, c)
% uses functions ring.m, projimg.m
 

% the names
names{1} = 'intensity';
names{2} = 'x-position';
names{3} = 'y-position';
names{4} = 'anisotropy';
names{5} = 'angle';
names{6} = 'volume outside cell';
names{7} = 'hole diameter';
names{8} = 'largest hole diameter';
names{9} = 'angle displacement';
names{10} = 'displacement';

% the units
units{1} = 'arbitrary';
units{2} = 'microns';
units{3} = 'microns';
units{4} = '';
units{5} = 'rad';
units{6} = 'microns^3';
units{7} = 'microns';
units{8} = 'microns';
units{9} = 'deg';
units{10} = 'microns';


%for myosin projection, only the top layer needs to be processed
% if c ~= 731 %358
% if c ~= 282; %1129
% % if (c <= 358 || c>369)
%      data = num2cell(NaN(10, 1));
%      return;
%  end

%keyboard

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
data{8} = NaN;
data{9} = NaN;
data{10} = NaN;





%% whole-cell properties

% get layers that were tracked; all others are NaN in x_values 
centroids = Cell.centroidStack(embryo.getCellStack(c, t)) * dx;
x_values = centroids(:, 2);
tracked = ~isnan(x_values);
ind_tracked = find(tracked);
ntracked = sum(tracked);



if (z ~= embryo.masterLayer || ntracked<20)
     return;
end


c
% get indicis of master, loest tracked, and highest layer
indml = embryo.translateZ(embryo.masterLayer)+1; % master layer
indlt = embryo.translateZ(embryo.lowestTracked(c,t))+1; % lowest tracked
indhl = length(x_values);                           % highest layer
nup = round(0.6/dz);     % select cell outline nup above lowest tracked
indnlt = ind_tracked(min([1+nup ntracked])); % index of nth lowest tracked

% get zs of several relevant layers
zlt = embryo.lowestTracked(c, t); % lowest tracked
znlt = embryo.unTranslateZ(indnlt-1);       % and nth lowest tracked

% for yolk, only the layers "nabv" above and
% "nbel" below the lowest layer will be processed
nabv=4;
nbel=6;
nabv = min([nabv indhl-indlt]);
nbel = min([nbel indlt-1]);

co_out = zeros(8,nabv+nbel+1);

%% get tilt at lower end of cells

% extract coordinates of centroids
y_values = centroids(:, 1);
x_values = centroids(:, 2);
z_values = (1:length(x_values)) * dz;
z_values = z_values(:);

% select coordinates of 'n_sel' closest layers to zlt that are not NaN
n_sel = min([round(5/dz) ntracked-nup]);    % number of z layers chosen to calculate tilt
ind_an = ind_tracked(1+nup:min([n_sel+nup ntracked]));



% compute tilt via principle components
x_sel = x_values(ind_an);
y_sel = y_values(ind_an);
z_sel = z_values(ind_an);
X = [x_sel,y_sel,z_sel];
[coeff,score,roots] = princomp(X);

% first component, vector pointing towards largest variation of centroids
v = coeff(:,1);
if (v(3)<0)   % and towards the top of the cell
    v=-v;
end

% fitted x and y centroid coordinates based on tilt
meanX = mean(X,1);
zc_fit = z_values(indlt-nbel:indlt+nabv);
z_offsets = meanX(1,3)-zc_fit;     % offset in microns
xc_fit=meanX(1,1)-z_offsets*v(1)/v(3);
yc_fit=meanX(1,2)-z_offsets*v(2)/v(3);

ang_tot = atan(sqrt(v(1)^2+v(2)^2)/abs(v(3)));
ang_dir = atan2(v(2),v(1));


%%

% first estimate background noise further away from yolk
yolkc_bg=[];
for ii=nabv+nbel+1+3:nabv+nbel+1+5
    ind_i = indlt-nbel-1+ii;
    z_i = embryo.unTranslateZ(ind_i-1);   % get index
    
    % get yolk for all cells
    yolk = getYolk(t, z_i);
%     if (tracked(ind_i))
%         [roic, R] = drawCellSmall(embryo, t, z_i, c);  % roic is cell_img
%         % (R(1) is offset in vertical, R(2) in horizontal direction)
%     else
        [roic, R] = drawCellSmall(embryo, t, znlt, c);  % roic is cell_img
        % (R(1) is offset in vertical, R(2) in horizontal direction)
        z_offset = z_values(indnlt)-z_values(ind_i);  % offset in microns
        xsh=round(-z_offset*v(1)/v(3)/dx);
        ysh=round(-z_offset*v(2)/v(3)/dx);
        R=R+[ysh xsh];
%    end
    yxs = size(roic);
    
    % extract yolk inside cell
    yolkc = double(yolk(R(1):R(1)+yxs(1)-1,R(2):R(2)+yxs(2)-1));
    yolkc = yolkc.*roic;
    yolkc = yolkc(:);
    yolkc_bg = [yolkc_bg yolkc'];
end

mnys = mean(yolkc_bg);
sdys = std(yolkc_bg);



for ii=1:nabv+nbel+1
    ind_i = indlt-nbel-1+ii;
    z_i = embryo.unTranslateZ(ind_i-1);   % get index
    
    % get yolk for all cells
    yolk = getYolk(t, z_i);
%     if (ind_i>indnlt && tracked(ind_i))
%         [roic, R] = drawCellSmall(embryo, t, z_i, c);  % roic is cell_img
%         % (R(1) is offset in vertical, R(2) in horizontal direction)
%     else
        [roic, R] = drawCellSmall(embryo, t, znlt, c);  % roic is cell_img
        % (R(1) is offset in vertical, R(2) in horizontal direction)
        z_offset = z_values(indnlt)-z_values(ind_i);  % offset in microns
        xsh=round(-z_offset*v(1)/v(3)/dx);      % xy shift in pixel
        ysh=round(-z_offset*v(2)/v(3)/dx);
        R=R+[ysh xsh];
%    end
    npad = 1;
    roic = padarray(roic,[npad npad]); % increase region around cell
    roic = bwmorph(roic,'thicken',npad);
    R = R-[npad npad];
    yxs = size(roic);
    
    
    % extract yolk inside cell
    yolkc = double(yolk(R(1):R(1)+yxs(1)-1,R(2):R(2)+yxs(2)-1));
    yolkc = yolkc.*roic;
    yolkc_raw = yolkc;
    yolkc_red = yolkc_raw - (mnys+sdys);
    yolkc_red(yolkc_red<0) = 0;
    intensity = mean(mean(yolkc_red(roic)));
    if (sum(sum(yolkc_red>0))<10)
       yolkc_red = yolkc_raw; 
       intensity = 0;
    end
    yolkc = yolkc_red/mean(mean(yolkc_red(roic))); % normalize
    %keyboard
    
    % project on plane perpendicular to cell axis
    roic = projimg(roic,ang_tot,ang_dir,...
        [xc_fit(ii)/dx-R(2) yc_fit(ii)/dx-R(1)]);
    yolkc = projimg(yolkc,ang_tot,ang_dir,...
        [xc_fit(ii)/dx-R(2) yc_fit(ii)/dx-R(1)]);
    yolkc = yolkc.*roic/cos(ang_tot);
    
    %keyboard
    
    %initial paramters for fitting:
    ext=(2^ceil(max(log2(size(yolkc)))));       %smooth yolkin inside cell
    [yolkf,filt] = get_filtered_gauss(yolkc,ext,1,0);
    [roicf,filt] = get_filtered_gauss(roic,ext,1,0);
    yolkf(roic>0)=yolkf(roic>0)./roicf(roic>0);
    yolkf=yolkf.*double(roic);
    
    yolk_th=yolkf>1;                    %and threshold
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
    
    %keyboard
    %     ng=length(mxax);
    
    xe0=0.2;      % set paramters and upper and lower bounds
    ub=1.0;
    lb=0;
    
    %xe0=[xe0 1 cents(1,1) cents(1,2) ...
    cenx = xc_fit(ii)/dx-R(2);
    ceny = yc_fit(ii)/dx-R(1);
    
    [indmx_y,indmx_x] = find(yolkf==max(yolkf(:)));
    angmx = atan2(indmx_y-ceny,indmx_x-cenx);
    
    maxradius=max([yxs(2) yxs(1)]);
    minx = max([cenx-2/dx 1]);
    maxx = min([cenx+2/dx yxs(2)]);
    miny = max([ceny-2/dx 1]);
    maxy = min([ceny+2/dx yxs(1)]);
    
    xe0=[xe0 4 cenx ceny max([maxradius/5 1/dx]) 1 orient 0.5/dx 0.5 angmx];
    ub=[ub 5 maxx maxy 0.4*maxradius 1.2  300 1.0/dx 1  300];
    lb=[lb 0 minx miny 0.0*maxradius 0.8 -300 0    0 -300];
    
    
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
        co_out(8,ii) = intensity;          % 'Stalk intensity'
    else
        co_out(:,ii)=NaN(8,1);
    end
    %keyboard
end

%% calculate properties that characterize yolk stalk

% pre selection 
min_ysint_fac = 0.66; %intensity factor below which images are discarded
[mxint,ind_mx_ys] = max(co_out(8,:));
ind_high_ys = find(co_out(8,:) > min_ysint_fac * max(co_out(8,:)));

% x and y position of YS bulk (intensity above min_ysint_fac) 
xyspot = [mean(co_out(1,ind_high_ys)) mean(co_out(2,ind_high_ys))];

data{2} = xyspot(1);
data{3} = xyspot(2);

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

xycent = [xc_fit(ind_mx_ys) yc_fit(ind_mx_ys)];
%[angle rel_dist] = calc_loc(roic,xycent,xyspot);

data{9}=atan2(xyspot(2)-xycent(2),xyspot(1)-xycent(1))/pi*180;
data{10}=norm(xyspot-xycent);

%keyboard
    



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

function [angle rel_dist] = calc_loc(roic,xycent,xyspot)
    % calculates location of myosin spot relative to cell centroid and boundary; 
    % outputs: orientation (angle) and relative distance from 
    % centroid (in [0 1])

    angle=atan2(xyspot(2)-xycent(2),xyspot(1)-xycent(1));

    a=size(roic);


    check=0;
    count=0;
    while (check==0)
        x=xycent+count*[cos(angle),sin(angle)];
        xr=round(x);
        if (xr(1)<1 | xr(1)>a(2) | xr(2)<1 | xr(2)>a(1))
            check=1;
        end
        if ~check
            if ~roic(xr(2),xr(1))
                check=1;
            else
                count=count+1;
            end
        end
    end

    rel_dist=sqrt(sum((xyspot-xycent).^2))/(count-1);
    rel_dist=min([rel_dist 1]);

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


