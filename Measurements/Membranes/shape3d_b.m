function [data names units] = shape3d_b(embryo, getMemb, t, z, c, dx, dz, dt, other)
% computes the tilt angles at (t, z) for Cell i given the Embryo4D,
% the membranes, and the relevant resolutions
% based on linear regresion in the vicinity of z

% the names
names{1} = 'Cell length';
names{2} = 'Volume';
names{3} = 'Centroids';
names{4} = 'Tilt vectors';
names{5} = 'All tilts-tot';
names{6} = 'All tilts-dir';
names{7} = 'Depth absolute';
names{8} = 'Depth relative';
names{9} = 'Volume accumulate';
names{10} = 'Cross section area';
names{11} = 'Surface_area';
names{12} = 'Touch_bot';
names{13} = 'Touch_top';
names{14} = 'Anisotropy';
names{15} = 'Angle anisotropy';
names{16} = 'Intensity all';
names{17} = 'Rough tilt vector';

% the units
units{1} = 'micron';
units{2} = 'micron^3';
units{3} = 'micron';
units{4} = 'micron';
units{5} = 'degree';
units{6} = 'degree';
units{7} = 'micron';
units{8} = '';
units{9} = 'micron^3';
units{10} = 'micron^2';
units{11} = 'micron^2';
units{12} = '';
units{13} = '';
units{14} = '';
units{15} = 'deg';
units{16} = '';
units{17} = '';
data = num2cell(NaN(length(units), 1));




%if c ~= 385; %358;% 340;%  676;% 209 %353
%if (c ~= 398 && c ~= 398 && c~=582  && c~=794)
% 815 756
%if (c~=582 && c~=978);%   537          815);  %784)
% if ( c ~= 385 && c~=784 && c ~= 794  && c~=815)
%if (c~=705 && c~= 442 && c~=582 && c~=200 && c~=815)
%  if (c ~= 385 && c~=784 && c ~= 794  && c~=815 && c ~= 789 && c~=682 && c ~= 575  && c~=443 && c~=1077)
if (t<=97 )%|| t~= 17)
   return;
end
% 


%% the following, do only for master layer since these properties are a 
% function of the whole cell stack, not an individual layer

if z ~= embryo.masterLayer
     return;
end

area = embryo.getCell(c, t, z).area * dx^2;
if (area>300)
   return; 
end

c


%% some definitions and indeces

% get layers that were tracked; all others are NaN in x_values 
centroids = Cell.centroidStack(embryo.getCellStack(c, t)) * dx;
x_values = centroids(:, 2);
y_values = centroids(:, 1);
tracked = ~isnan(x_values);
ind_tracked = find(tracked);
ntracked = sum(tracked);

% get indicis of master, highest tracked, and highest layer
indml = embryo.translateZ(embryo.masterLayer)+1; % master layer
indht = embryo.translateZ(embryo.highestTracked(c,t))+1; % highest tracked
indlt = embryo.translateZ(embryo.lowestTracked(c,t))+1; % highest tracked
indhl = length(x_values);                           % highest layer
indll = 1;
indi = embryo.translateZ(z)+1;  % index of current layer





%% estimate tilts and centroids for all layers

if isnan(dt)
    % number of z layers chosen to calculate tilt
    dz_sel = 6;
    n_sel = round(dz_sel/dz);
    
    % avoid layers at bottom and top
    dz_avoid = 3;
    n_avoid = round(dz_avoid/dz);
else 
    % number of z layers chosen to calculate tilt
    dz_sel = 5;
    n_sel = ceil(dz_sel/dz);
    
    % avoid layers at bottom and top
    dz_avoid = 2;
    n_avoid = round(dz_avoid/dz);
    t
end    
    

if (ntracked<(n_sel+1+2*n_avoid));  % exit if not enough layers available
    return;
end

ind_tracked_avoid = ind_tracked(n_avoid+1:end-n_avoid);
indlb = ind_tracked_avoid(1);   % lower and 
indub = ind_tracked_avoid(end);  % lower end of bulk of cell

centroids = Cell.centroidStack(embryo.getCellStack(c, t)) * dx;
y_values = centroids(:, 1);
x_values = centroids(:, 2);
z_values = (1:length(x_values)) * dz;
z_values = z_values(:);

%keyboard

Xests = zeros(indhl,3);
tilt_tots = zeros(indhl,1);
tilt_dirs = zeros(indhl,1);
vs = zeros(indhl,3);
for i=1: indhl
    indi = i;
    % select coordinates of closest layers to z that are not NaN
    [q,ind_close] = sort(abs(indi-ind_tracked_avoid));
    ind_an = ind_tracked_avoid(ind_close(1:n_sel));
    
    % compute tilt via principle components
    X = [x_values,y_values,z_values];
    X = X(ind_an,:);
    [coeff,score,roots] = princomp(X);
    % first component, vector pointing towards largest variation of centroids
    v = coeff(:,1);
    if (v(3)<0)   % and towards the top of the cell
        v=-v;
    end
    
    % tilt angles
    tilt_tots(i) = atan(sqrt(v(1)^2+v(2)^2)/abs(v(3)));
    tilt_dirs(i) = atan2(v(2),v(1));
    
    % estimate centroid
    Xm = mean(X,1);  % mean of centroids used for fitting
    zsh = z_values(indi)-Xm(3);  % offset to current depth, in microns
    xsh=zsh*v(1)/v(3);   % x and y shift in micron
    ysh=zsh*v(2)/v(3);
    Xests(i,1) = Xm(1)+xsh;
    Xests(i,2) = Xm(2)+ysh;
    Xests(i,3) = z_values(indi);
    vs(i,:) = v;
end

%% test run to get intensity (for live data)

if ~isnan(dt)
    intens_alll = zeros(indhl,1);
    
    intens_diffs=[];
    
    roi0 = 0*getMemb(t, embryo.masterLayer);  % zero background array
    siz = size(roi0);
    
    i=1;
    while  (i<=indhl)
        z_i = embryo.unTranslateZ(i-1);  % index of layer
        intens = getMemb(t, z_i);  % global memb at layer zub
        [roicn Rn] = get_roi_combined(embryo,roi0,i,ind_tracked,4,2,vs(i,:),dx,dz,c,t);
        yxsn = size(roicn);
        % split roi up into interior and exterior
        roic_centern = bwmorph(roicn,'thin',ceil(min(yxsn)/2/3));
        roic_bordern = roicn-roic_centern;
        %
        intensc_icent = sum(sum(roic_centern.*intens(Rn(1):Rn(1)+yxsn(1)-1,Rn(2):Rn(2)+yxsn(2)-1)));
        intensc_ibord = sum(sum(roic_bordern.*intens(Rn(1):Rn(1)+yxsn(1)-1,Rn(2):Rn(2)+yxsn(2)-1)));
        intens_diff = intensc_ibord/sum(roic_bordern(:))-...
            intensc_icent/sum(roic_centern(:));
        intens_diffs = [intens_diffs intens_diff];
        %intens_all(i) = intensc_ibord/sum(roic_bordern(:));
        intens_alll(i) = intens_diff;
        i=i+1;
    end
    [mxalll,indmxalll] = max(intens_alll);
    minalll1 = max([0 min(intens_alll(1:indmxalll))]);
    minalll2 = max([0 min(intens_alll(indmxalll:end))]);
end





%% upper end of cell: estimate length, volume, and maximal centroid

if isnan(dt)
    th_diff = 2;
    i=indub-2;
else
    th_diff = minalll2 + (mxalll-minalll2)/6;  %15
    i=indub-4;
end

intens_all = zeros(indhl,1); 

intens_diffs=[];
intens_diff=th_diff+10;

roi0 = 0*getMemb(t, embryo.masterLayer);  % zero background array
siz = size(roi0);

%i=1;
while  (i<=indhl && intens_diff>=th_diff)  
    z_i = embryo.unTranslateZ(i-1);  % index of layer
    intens = getMemb(t, z_i);  % global memb at layer zub
    [roicn Rn] = get_roi_combined(embryo,roi0,i,ind_tracked,4,2,vs(i,:),dx,dz,c,t);
    yxsn = size(roicn);
    % split roi up into interior and exterior
    roic_centern = bwmorph(roicn,'thin',ceil(min(yxsn)/2/3));
    roic_bordern = roicn-roic_centern;
    % 
    intensc_icent = sum(sum(roic_centern.*intens(Rn(1):Rn(1)+yxsn(1)-1,Rn(2):Rn(2)+yxsn(2)-1)));
    intensc_ibord = sum(sum(roic_bordern.*intens(Rn(1):Rn(1)+yxsn(1)-1,Rn(2):Rn(2)+yxsn(2)-1)));
    intens_diff = intensc_ibord/sum(roic_bordern(:))-...
                  intensc_icent/sum(roic_centern(:));
    intens_diffs = [intens_diffs intens_diff];    
    %intens_all(i) = intensc_ibord/sum(roic_bordern(:));
    intens_all(i) = intens_diff;
    i=i+1;
end
i=i-2;

intens_up = intens_diffs;
%keyboard

if (i<indht)
    centroids(i+1:end,:) = NaN;
    x_values = centroids(:, 2);
    y_values = centroids(:, 1);
    tracked = ~isnan(x_values);
    ind_tracked = find(tracked);
    ntracked = sum(tracked);
    % get indicis of master, highest tracked, and highest layer
    if ntracked>0    
        indht = ind_tracked(end); % highest tracked
    else
        return;
    end
    
end

% monitor whether cell touches uper layer 
touch_top = (indht==indhl);

% if ~isnan(dt)
%     i=min([i+1 indhl]);
% end
Xub = Xests(i,:);


if (ntracked<(n_sel+2*n_avoid));  % exit if not enough layers available
    return;
end


%% lower end of cell: estimate length, volume, and maximal centroid 

if isnan(dt)
    th_diff = 2;
    i=indlb+3;
elseif dz==1
    th_diff = minalll1+(mxalll-minalll1)/4;  % 20;
    i=indlb+4;
else
    th_diff = 0;  % 20;   
    i=indlb+4;
end


intens_diffs=[];
intens_diff=th_diff+10;
intens_diff_pre=[th_diff+10];

% extrapolate
%while  (i>=indll && (abs(intens_diff)>=2 || abs(intens_diff_pre)>=2))
while  (i>=indll && intens_diff>=th_diff)% || intens_diff_pre>=2)
    intens_diff_pre = intens_diff;
    z_i = embryo.unTranslateZ(i-1);  % index of layer
    intens = getMemb(t, z_i);  % global memb at layer zlb
    [roicn Rn] = get_roi_combined(embryo,roi0,i,ind_tracked,4,2,vs(i,:),dx,dz,c,t);
    yxsn = size(roicn);
    % split roi up into interior and exterior
    roic_centern = bwmorph(roicn,'thin',ceil(min(yxsn)/2/3));
    roic_bordern = roicn-roic_centern;
    % 
    intensc_icent = sum(sum(roic_centern.*intens(Rn(1):Rn(1)+yxsn(1)-1,Rn(2):Rn(2)+yxsn(2)-1)));
    intensc_ibord = sum(sum(roic_bordern.*intens(Rn(1):Rn(1)+yxsn(1)-1,Rn(2):Rn(2)+yxsn(2)-1)));
    intens_diff = intensc_ibord/sum(roic_bordern(:))-...
                  intensc_icent/sum(roic_centern(:));
    intens_diffs = [intens_diffs intens_diff];    
    intens_all(i) = intens_diff;

    i=i-1;
    %keyboard
end
i=i+2;

%keyboard

% monitor whether cell touches lowest layer 
touch_bot = (indlt==indll);

% another one up
%i=i+1;


intens_dow = intens_diffs;

if (i>indlt)
    centroids(1:i-1,:) = NaN;
    x_values = centroids(:, 2);
    y_values = centroids(:, 1);
    tracked = ~isnan(x_values);
    ind_tracked = find(tracked);
    ntracked = sum(tracked);
    % get indicis of master, highest tracked, and highest layer
    if ntracked>0    
        indlt = ind_tracked(1); % highest tracked
    else
        return;
    end        
end

% if ~isnan(dt)
%     i=max([i-1 indll]);
% end
Xlb = Xests(i,:);







%keyboard



%% repair cells
%% some definations

navoid = 0 ;            % avoid n layers at end
indhta = ind_tracked(end-navoid);
indlta = ind_tracked(1+navoid);

dz_extrap = 10;   % distance of extrapolation (in micron)
nextrap = round(dz_extrap/dz);

dz_back = 1.8;   % distance of points to look back, for extrapolation (in micron)
nback = round(dz_back/dz);

nclose = max([3 ceil(3.6/dz)]);  % number of layers to inter and extrapolate cents, verts, and areas (for cell reconstruction)

%nroi = ceil(0.9/dz);  % number of rois to combine at ends for extrapolation


if (ntracked<(nclose+2));  % exit if not enough layers available
    return;
end


%% collect, interpolate and extrapolate all centroids and vertices  

cents = NaN(indhl,3);  % default centroids; important if master layer at top end
areas = NaN(indhl,1);
verts = cell(indhl,1);
col = NaN(indhl,1);

% march from bottom to top end and insert existing values
for i=indll:indhl
    zind = embryo.unTranslateZ(i-1);   % index of layer
    if tracked(i)   
        % extract coordinates of centroids
        areas(i) = embryo.getCell(c, t, zind).area * dx^2;
        cents(i,:) = [x_values(i) y_values(i) i*dz];
        vert = fliplr(embryo.getCell(c, t, zind).vertexCoords * dx);
        verts{i} =  [vert i*dz+zeros(size(vert,1),1)];
        col(i) = 1;
    end
end

%keyboard
% repair nonexisting values
for i=indll:indhl
    if ~tracked(i)   
        [q,ind_sort] = sort(abs(i-ind_tracked));
        ind_close = ind_tracked(ind_sort(1:nclose));
        pfit_cx = polyfit(ind_close,cents(ind_close,1),1);
        pfit_cy = polyfit(ind_close,cents(ind_close,2),1);
        pfit_ar = polyfit(ind_close,areas(ind_close),1);
        
        cents(i,1) = polyval(pfit_cx,i);
        cents(i,2) = polyval(pfit_cy,i);
        cents(i,3) = i*dz;
        areas(i) = polyval(pfit_ar,i);

        indc = ind_close(1);
        vert = verts{indc};
        arfac = areas(indc)/areas(i);
        arfac = max([0.2 arfac]);
        vts = [vert(:,1)-cents(indc,1) vert(:,2)-cents(indc,2)];
        vts = [vts(:,1)/sqrt(arfac)...
               vts(:,2)/sqrt(arfac)];
        verts{i} = [vts(:,1)+cents(i,1) vts(:,2)+cents(i,2) vts(:,2)*0+i*dz];
%        keyboard
        col(i) = 0.5;        
    end
end 

vertices = [];
vertices_sel = [];
for i=indll:indhl
    vertices = [vertices;verts{i}];
    if (i>=indlta && i<=indhta)
       vertices_sel = [vertices_sel;verts{i}]; 
    end
end




% interpolate vertices for interpolating cell
dzz = 2*dx; %new resolution in z 

zzs = dzz:dzz:(indhl*dz+dzz*0.99);
nzzs = length(zzs);
ind_dzz = find(zzs>indlta*dz & zzs<indhta*dz);
indlta_dzz = ind_dzz(1)-1;
indhta_dzz = ind_dzz(end)+1;

verts_dzz = cell(nzzs,1);

vertices_dzz = [];
for ii=indlta_dzz:indhta_dzz;
    zz=zzs(ii);
    [q,ind_sort] = sort(abs(zz/dz-ind_tracked));
    ind_close = ind_tracked(ind_sort(1:nclose));
    pfit_cx = polyfit(ind_close*dz,cents(ind_close,1),1);
    pfit_cy = polyfit(ind_close*dz,cents(ind_close,2),1);
    pfit_ar = polyfit(ind_close*dz,areas(ind_close),1);
    
    cents_x = polyval(pfit_cx,zz);
    cents_y = polyval(pfit_cy,zz);
    areas_dzz = polyval(pfit_ar,zz);
    
    indc = ind_close(1);
    vert = verts{indc};
    arfac = areas_dzz/areas(indc);
    arfac = max([0.2 arfac]);
    vts = [vert(:,1)-cents(indc,1) vert(:,2)-cents(indc,2)];
    vts = [vts(:,1)/sqrt(arfac)...
        vts(:,2)/sqrt(arfac)];
    vert_dzz = [vts(:,1)+cents_x vts(:,2)+cents_y vts(:,2)*0+zz];
    %        keyboard
    vertices_dzz = [vertices_dzz;vert_dzz];
    verts_dzz{ii} = vert_dzz;
end 






%keyboard







%% repair upper and lower end; assume missing parts are convex


%upper end
if isnan(dt)
    memth = 0.5;
else
    memth = 0.0;  % threshold for elimation of points
end


% estimate area and anisotpy of rois (for extrapolation)
majors = [];
minors = [];
orients = [];
zs = [];
ars = [];
for i=indhta-nback+1:indhta
    %keyboard
    if isnan(dt)
        v = vs(i,:);
    else
        v = get_tiltv(cents,i,nclose);  % get tilt vector
    end
    %keyboard
    [roic R] = get_roi_combined(embryo,roi0,i,ind_tracked,4,2,v,dx,dz,c,t); % get combind roi
    z_i = embryo.unTranslateZ(i-1);
    props = regionprops(roic,'Area','Centroid','MajorAxisLength', 'MinorAxisLength', 'Orientation');
    majors = [majors props(1).MajorAxisLength*dx];
    minors = [minors props(1).MinorAxisLength*dx];
    orients = [orients -props(1).Orientation/180*pi];
    zs = [zs i*dz];
    ars = [ars props(1).Area*dx^2];
end
centroids = [props.Centroid];
nspots = size(centroids,2)/2;
centroids = transpose(reshape(centroids,[2 nspots]));


%fitting properties 
orient = angle(mean(exp(complex(0,1)*orients)));
if (length(majors)==1)
    majors = [majors majors];
    minors = [minors minors];
    ars = [ars ars];
    zs = [zs-dz zs];
end
pfit_mj = polyfit(zs,majors,1);
pfit_mn = polyfit(zs,minors,1);
pfit_ars = polyfit(zs,ars,1);


%keyboard

% collect points (of high membrane intensity) across top layers
yxs = size(roic);
Rm=R;
vertices_up = [];
membcc = roic;
i=indhta;
while ( (i<=min([indhta+nextrap indhl]) && (sum(sum(membcc))>0) && ...
         polyval(pfit_ars,i*dz)>0) || i==indhta)
        %polyval(pfit_mj,i*dz)>0 && polyval(pfit_mn,i*dz)>0)
    %keyboard
    z_offset = (i-indhta)*dz;  % offset in microns
    xsh=round(z_offset*v(1)/v(3)/dx);      % xy shift in pixel
    ysh=round(z_offset*v(2)/v(3)/dx);
    R=Rm+[ysh xsh];
    

    zind = embryo.unTranslateZ(i-1);   % index of layer
    memb = getMemb(t, zind);  % global memb at layer zub
    %memb = other.YolkStalk(t, zind);
  
    R = [min([R(1) siz(1)-yxs(1)+1]) min([R(2) siz(2)-yxs(2)+1])];
    R = [max([R(1) 1]) max([R(2) 1])];

    membc = double(memb(R(1):R(1)+yxs(1)-1,R(2):R(2)+yxs(2)-1));
    
    ext=(2^ceil(max(log2(size(membc)))));       %smooth yolkin inside cell
    [membf,filt] = get_filtered_gauss(membc,ext,1,0);
    [membf_norm,filt] = get_filtered_gauss(membc*0+1,ext,1,0);
    membf=membf./membf_norm;

    % resize roi
    %mj = polyval(pfit_mj,i*dz)/polyval(pfit_mj,indhta*dz);
    %mn = polyval(pfit_mn,i*dz)/polyval(pfit_mn,indhta*dz);
    %roicc = scaleimg(roic,mj,mn,orient,centroids);
    
    ar = polyval(pfit_ars,i*dz)/polyval(pfit_ars,indhta*dz);
    roicc = scaleimg(roic,sqrt(ar),sqrt(ar),orient,centroids);
    

    
    
    if ~isnan(dt)
        % split roi up into interior and exterior
        yxsn = size(roicc);
        roic_centern = bwmorph(roicc,'thin',ceil(min(yxsn)/2/2));
        %roic_centern = roicc;
        %roic_centern = bwmorph(roic_centern,'skel',1);
        roic_bordern = roicc-roic_centern;
        membff = roic_centern.*membf;
        intensc_icent = mean(membff(roic_centern==1))+...
                           0*std(membff(roic_centern==1));
        membf=membf-intensc_icent;
    end

    membcc = membf.*roicc;
    membcc(membcc<0)=0;
   
    
    if (i==indhta)
        roic_up=roic;
        mxm=max(membcc(:));
        [yk,xk] = find(roic);
    else
        membcc = membcc>memth*mxm;
        if ~isnan(dt)
            props = regionprops(membcc,'Area');
            ars = [props.Area]*dx^2;
            %if sum(membcc(:)) < 3/dx^2
            if max(ars) < 1.7
       %         membcc=membcc*0;
            end
        end
        [yk,xk] = find(membcc);        
    end
    
    if ~isempty(xk)
       points = [dx*(xk+R(2)),dx*(yk+R(1)),dz*i*ones(size(xk,1),1)];
       vertices_up = [vertices_up;points];
    end
    %keyboard
    i=i+1;
end

%proj_up = (repmat(Xub,size(vertices_up,1),1)-vertices_up)*transpose(vs(indhta,:));
%vertices_up(proj_up<0,:) = [];
ind_toohigh = vertices_up(:,3)>Xub(3);
vertices_up(ind_toohigh,:) = [];


%keyboard




%lower end
if isnan(dt)
    memth = 0.7;
else
    memth = 0.0;
end

% estimate anisotpy of rois (for extrapolation)
majors = [];
minors = [];
orients = [];
ars= [];
zs = [];
for i=indlta+nback-1:-1:indlta
    %v = get_tiltv(cents,i,nclose);  % get tilt vector
    if isnan(dt)
        v = vs(i,:);
    else
        v = get_tiltv(cents,i,nclose);  % get tilt vector
    end
    [roic R] = get_roi_combined(embryo,roi0,i,ind_tracked,4,2,v,dx,dz,c,t); % get combind roi
    z_i = embryo.unTranslateZ(i-1);
    props = regionprops(roic,'Area','Centroid','MajorAxisLength', 'MinorAxisLength', 'Orientation');
    majors = [majors props(1).MajorAxisLength*dx];
    minors = [minors props(1).MinorAxisLength*dx];
    orients = [orients -props(1).Orientation/180*pi];
    ars = [ars props(1).Area*dx^2];
    zs = [zs i*dz];
end
centroids = [props.Centroid];
nspots = size(centroids,2)/2;
centroids = transpose(reshape(centroids,[2 nspots]));


%fitting properties 
orient = angle(mean(exp(complex(0,1)*orients)));
if (length(majors)==1)
    majors = [majors majors];
    minors = [minors minors];
    ars = [ars ars];
    zs = [zs-dz zs];
end
pfit_mj = polyfit(zs,majors,1);
pfit_mn = polyfit(zs,minors,1);
pfit_ars = polyfit(zs,ars,1);




% collect points (of high membrane intensity) across bottom layers
yxs = size(roic);
vertices_dow = [];
Rm=R;
membcc = roic;
i=indlta;
while ( (i>=max([indlta-nextrap indll]) && (sum(sum(membcc))>0) && ...
        polyval(pfit_ars,i*dz)>0) || i==indlta)
        %   polyval(pfit_mj,i*dz)>0 && polyval(pfit_mn,i*dz)>0)
    %keyboard
    z_offset = (i-indlta)*dz;  % offset in microns
    xsh=round(z_offset*v(1)/v(3)/dx);      % xy shift in pixel
    ysh=round(z_offset*v(2)/v(3)/dx);
    R=Rm+[ysh xsh];
    zind = embryo.unTranslateZ(i-1);   % index of layer
    memb = getMemb(t, zind);  % global memb at layer zub
    %memb = other.YolkStalk(t, zind);
    R = [min([R(1) siz(1)-yxs(1)+1]) min([R(2) siz(2)-yxs(2)+1])];
    R = [max([R(1) 1]) max([R(2) 1])];

    membc = double(memb(R(1):R(1)+yxs(1)-1,R(2):R(2)+yxs(2)-1));    
    ext=(2^ceil(max(log2(size(membc)))));       %smooth intensity inside cell
    [membf,filt] = get_filtered_gauss(membc,ext,1,0);
    [membf_norm,filt] = get_filtered_gauss(membc*0+1,ext,1,0);
    membf=membf./membf_norm;

    % resize roi
%     mj = polyval(pfit_mj,i*dz)/polyval(pfit_mj,indlta*dz);
%     mn = polyval(pfit_mn,i*dz)/polyval(pfit_mn,indlta*dz);
%     roicc = scaleimg(roic,mj,mn,orient,centroids);
    ar = polyval(pfit_ars,i*dz)/polyval(pfit_ars,indlta*dz);
    roicc = scaleimg(roic,sqrt(ar),sqrt(ar),orient,centroids);

%    keyboard
%     if ~isnan(dt)
%         membf=membf-mean(membf(:));
%     end    

    if ~isnan(dt)
        % split roi up into interior and exterior
        yxsn = size(roicc);
        %roic_centern = bwmorph(roicc,'shrink',ceil(min(yxsn)/2/4));
        roic_centern = bwmorph(roicc,'thin',ceil(min(yxsn)/2/2));
        roic_bordern = roicc-roic_centern;
        membff = roic_centern.*membf;
        intensc_icent = mean(membff(roic_centern==1))+...
                           0*std(membff(roic_centern==1));
        membf=membf-intensc_icent;
    end

    membcc = membf.*roicc;
    membcc(membcc<0)=0;


    if (i==indlta)
        roic_dow=roic;
        mxm=max(membcc(:));
        [yk,xk] = find(roic);
    else
        membcc = membcc>memth*mxm;
        if ~isnan(dt)
            props = regionprops(membcc,'Area');
            ars = [props.Area]*dx^2;
            %if sum(membcc(:)) < 3/dx^2
            if max(ars) < 2
          %      membcc=membcc*0;
            end
        end
        [yk,xk] = find(membcc);        
    end
    
    %keyboard
    if ~isempty(xk)
       points = [dx*(xk+R(2)),dx*(yk+R(1)),dz*i*ones(size(xk,1),1)];
       vertices_dow = [points;vertices_dow];
    end
    i = i-1;
end
%keyboard
% eliminate points beyond intensity TH
% proj_dow = (repmat(Xlb,size(vertices_dow,1),1)-vertices_dow)*transpose(vs(indlta,:));
% vertices_dow(proj_dow>0,:) = [];
ind_toolow = vertices_dow(:,3)<Xlb(3);
vertices_dow(ind_toolow,:) = [];







%% volume bulk

% size of region
vertices_all = [vertices_dzz;vertices_dow;vertices_up];

% limit values of cell 
lims = volumebounds(vertices_all(:,1),vertices_all(:,2),vertices_all(:,3),ones(size(vertices_all,1),1));
lims = [lims(1)-mod(lims(1),dx),lims(2)-mod(lims(2),dx)+dx, lims(3)-mod(lims(3),dx),lims(4)-mod(lims(4),dx)+dx, lims(5),lims(6)];

                
[qx,qy] = meshgrid(lims(1):dx:lims(2),lims(3):dx:lims(4));                  
               
vol = zeros(size(qx,1),size(qx,2),length(zzs));

for i=indlta_dzz:indhta_dzz
    vert = verts_dzz{i};
    vol(:,:,i) = inpolygon(qx,qy,vert(:,1),vert(:,2));
end

vol_bulk = vol;


                  
%keyboard
                  
 
                  

%% interpolate uper end 

% upper part
% range of data
lims_up = volumebounds(vertices_up(:,1),vertices_up(:,2),vertices_up(:,3),ones(size(vertices_up,1),1));
lims_up = [lims_up(1)-mod(lims_up(1),dx),lims_up(2)-mod(lims_up(2),dx)+dx, lims_up(3)-mod(lims_up(3),dx),lims_up(4)-mod(lims_up(4),dx)+dx, lims_up(5),lims_up(6)];
                
% arays containing coordinates                
indhl_dzz = find(zzs>lims_up(6));
indhl_dzz = indhl_dzz(1);
[qx_up,qy_up,qz_up] = meshgrid(lims_up(1):dx:lims_up(2),lims_up(3):dx:lims_up(4),...
                          zzs(indhta_dzz):dzz:zzs(indhl_dzz));                  
sz_up = size(qx_up);
q_up = [qx_up(:) qy_up(:) qz_up(:)];

%keyboard
% generate convex hull of vertices
if single(max(vertices_up(:,3)))==single(indht*dz)
    FF = TriScatteredInterp(vertices_up(:,1:2),ones(size(vertices_up,1),1));
    layer=FF(qx_up,qy_up);
    vol_up = ~isnan(layer);
    sz_up = [sz_up 1];
else
    K_up = convhulln(vertices_up);
    in_up = inhull(q_up,vertices_up,K_up);
    vol_up = reshape(in_up,sz_up(1),sz_up(2),sz_up(3));
end


% lower part
% range of data
lims_dow = volumebounds(vertices_dow(:,1),vertices_dow(:,2),vertices_dow(:,3),ones(size(vertices_dow,1),1));
lims_dow = [lims_dow(1)-mod(lims_dow(1),dx),lims_dow(2)-mod(lims_dow(2),dx)+dx, lims_dow(3)-mod(lims_dow(3),dx),lims_dow(4)-mod(lims_dow(4),dx)+dx, lims_dow(5),lims_dow(6)];

% arays containing coordinates                
indll_dzz = find(zzs<lims_dow(5));
indll_dzz = indll_dzz(end);
[qx_dow,qy_dow,qz_dow] = meshgrid(lims_dow(1):dx:lims_dow(2),lims_dow(3):dx:lims_dow(4),...
                                  zzs(indll_dzz):dzz:zzs(indlta_dzz));   
       
sz_dow = size(qx_dow);
q_dow = [qx_dow(:) qy_dow(:) qz_dow(:)];

% generate convex hull of vertices
if single(min(vertices_dow(:,3)))==single(indlt*dz)
    FF = TriScatteredInterp(vertices_dow(:,1:2),ones(size(vertices_dow,1),1));
    layer=FF(qx_dow,qy_dow);
    vol_dow = ~isnan(layer);
    sz_dow = [sz_dow 1];    
else
    K_dow = convhulln(vertices_dow);
    in_dow = inhull(q_dow,vertices_dow,K_dow);
    vol_dow = reshape(in_dow,sz_dow(1),sz_dow(2),sz_dow(3));
end


%keyboard




%% put cell together 

% index of first element
Rx = find(single(lims(1):dx:lims(2))==single(qx_up(1)));
Ry = find(single(lims(3):dx:lims(4))==single(qy_up(1))); 
R3_up = [max([Ry 2]),max([Rx 2]),indhta_dzz];

Rx = find(single(lims(1):dx:lims(2))==single(qx_dow(1)));
Ry = find(single(lims(3):dx:lims(4))==single(qy_dow(1)));
R3_dow = [max([Ry 2]),max([Rx 2]),indlta_dzz];

vol_tot = vol_bulk;
vol_tot(:,:,[indlta_dzz indhta_dzz]) = 0;

%keyboard
vol_tot(R3_up(1)-1:R3_up(1)+sz_up(1)-2,R3_up(2)-1:R3_up(2)+sz_up(2)-2,R3_up(3):R3_up(3)+sz_up(3)-1) = vol_up;
vol_tot(R3_dow(1)-1:R3_dow(1)+sz_dow(1)-2,R3_dow(2)-1:R3_dow(2)+sz_dow(2)-2,R3_dow(3)-sz_dow(3)+1:R3_dow(3)) = vol_dow;



                  
%keyboard
clear vol_up vol_dow vol_bulk; %vol intens memb roi0;

%keyboard
% if (t==103 && c== 61)
%    keyboard;
% end



%%  calculat midline of cell

slw = 2;  % width of slice
nfit = max([3 ceil(3/dz)]);  % number of cents to fit tilt vector

% coordinates for entire cell
%[qx,qy,qz] = meshgrid(lims(1):dx:lims(2),lims(3):dx:lims(4),lims(5):dz:lims(6));                  
[qx,qy,qz] = meshgrid(lims(1):dx:lims(2),lims(3):dx:lims(4),zzs);                  

% points in cell 
rcell = find(vol_tot); 
Xcell = [qx(rcell) qy(rcell) qz(rcell)];
nrcell = length(rcell);

%keyboard
nitadj=1;
XXf = Xests;
VV = vs;
for i=1:nitadj  
    [XX,XXf,VV] = adjust_center_line(XXf,VV,Xcell,slw,nfit,dz);
end








%% display everything


plott = 0;
printyes=0;  % save plots


if plott

% centroids and vertices
figure(1);clf;
plot3(vertices(:,1),vertices(:,2),vertices(:,3),'r*');
hold on;
scatter3(cents(:,1),cents(:,2),cents(:,3),20,col);
hold off;
axis equal;

figure(2);clf;
plot3(vertices_sel(:,1),vertices_sel(:,2),vertices_sel(:,3),'r*');

figure(12);clf;
scatter3(cents(:,1),cents(:,2),cents(:,3),20,col);
axis equal;

% plot points of upper and lower end
figure(3);clf;
plot3(vertices_up(:,1),vertices_up(:,2),vertices_up(:,3),'o','LineWidth',1);
axis equal;
figure(4);clf;
plot3(vertices_dow(:,1),vertices_dow(:,2),vertices_dow(:,3),'o');
axis equal;

% display bulk of cell                  
fv = isosurface(qx,qy,qz,vol_bulk,0.5);
figure(7);clf;
pch = patch(fv,'FaceColor','green','EdgeColor','none');
%isonormals(qx,qy,qz,vol,pch)
view(3); daspect([1,1,1]); 
axis equal
camlight; camlight(-80,-10); lighting flat;

%keyboard
% display both parts

if (sz_up(3)~=1)
    fv = isosurface(qx_up,qy_up,qz_up,vol_up,0.5);
    figure(10);clf;
    pch = patch(fv,'FaceColor','green','EdgeColor','non');
    axis equal
    camlight; camlight(-80,-10); lighting flat;
    lighting gouraud
end

if (sz_dow(3)~=1)
    fv = isosurface(qx_dow,qy_dow,qz_dow,vol_dow,0.5);
    figure(11);clf;
    pch = patch(fv,'FaceColor','green','EdgeColor','non');
    axis equal
    camlight; camlight(-80,-10); lighting flat;
    lighting gouraud
end

% whole cell
fv = isosurface(qx,qy,qz,vol_tot,0.5);
figure(13);clf;
pch = patch(fv,'FaceColor','green','EdgeColor','none');
isonormals(qx,qy,qz,vol_tot,pch)
%view(3); daspect([1,1,1]); 
axis equal
camlight; camlight(-80,-10); lighting flat;
extr = [Xub;Xlb];
hold on;
plot3(extr(:,1),extr(:,2),extr(:,3),'ko','LineWidth',2,'MarkerSize',20);%,'MarkerFaceColor','k')
hold off;
axis([lims(1:4) Xlb(3)-3 Xub(3)+3])



figure(19);clf;
fv = isosurface(qx,qy,qz,vol_tot,0.5);
pch = patch(fv,'FaceColor','green','EdgeColor','none');
isonormals(qx,qy,qz,vol_tot,pch)
axis equal
camlight; camlight(-80,-10); lighting flat;

indXX = find(~isnan(XX(:,1))); 
extr = [XX(indXX(1),:);XX(indXX(end),:)];
hold on;
plot3(extr(:,1),extr(:,2),extr(:,3),'ko','LineWidth',2,'MarkerSize',20);%,'MarkerFaceColor','k')
plot3(XX(:,1),XX(:,2),XX(:,3),'k','LineWidth',2);%,'MarkerFaceColor','k')
hold off;
%axis([lims(1:4) Xlb(3)-3 Xub(3)+3])

figure(20);clf;
plot3(XXf(:,1),XXf(:,2),XXf(:,3),'k','LineWidth',2);%,'MarkerFaceColor','k')
axis equal

figure(21);clf;
plot3(XX(:,1),XX(:,2),XX(:,3),'k','LineWidth',2);%,'MarkerFaceColor','k')
axis equal

hold on;
plot3(Xests(:,1),Xests(:,2),Xests(:,3),'g','LineWidth',2);%,'MarkerFaceColor','k')
hold off;
axis equal
axis([lims(1:4) Xlb(3)-3 Xub(3)+3])


hh = figure(22);clf;
fv = isosurface(qx,qy,qz,vol_tot,0.5);
pch = patch(fv,'FaceColor','green','EdgeColor','none');
isonormals(qx,qy,qz,vol_tot,pch)
view(0,10);
axis equal
camlight; camlight(-80,-10); lighting flat;
%keyboard
if (printyes==1)
    print(hh, '-depsc2',['../../../../4d_analysis/live/cells_2ph-cad2/c' int2str(c) '_t' int2str(t) '.eps']);
end

end


%% several useful index scalars and vectors

%index of basal and apical end of midline of cell
indcell = find(~isnan(XX(:,1)));
indbe = indcell(1);
indae = indcell(end);

%keyboard
%% calculate quantities

indcell = find(~isnan(XX(:,1)));

% volume
cell_volume = nrcell*dx^2*dzz;

% length
cell_length = sum(sqrt(sum(diff(XX(indcell,:)).^2,2)));

% cross section area and cumulative volume
slw = 2;
[csareas,cumvols,anisos,ang_anisos] = calc_csarea_cumvol(XX,VV,Xcell,slw,dzz,dx,indbe,indae);


% depths
dists = sqrt(sum(diff(XX,1).^2,2));
dists(isnan(dists)) = 0;
dists = [dists;0];

% cumulative sum of distances
depths_abs = flipud(cumsum(flipud(dists)));
depths_abs = depths_abs-depths_abs(indae);
depths_abs(depths_abs<0) = 0;
%depths_rel = depths_abs/max(depths_abs);

% distances extended to lowest and highest layer
depths_ext = dists;
depths_ext(indll:indbe-1)=dists(indbe);
depths_ext(indae:indhl)=dists(indae-1);
depths_abs_ext = flipud(cumsum(flipud(depths_ext)));
depths_abs_ext = depths_abs_ext-depths_abs_ext(indae);
depths_rel_ext = depths_abs_ext/max(depths_abs);

tilt_tots = zeros(indhl,1);
tilt_dirs = zeros(indhl,1);
for i=1:indhl
    v = VV(i,:);
    tilt_tots(i) = atan(sqrt(v(1)^2+v(2)^2)/abs(v(3)));
    tilt_dirs(i) = atan2(v(2),v(1));
end

%% surface of cell

% side 
indlayers = [];
perims = zeros(indhl,1);
for i=1:nzzs
   slice2d = vol_tot(:,:,i);
   if (sum(slice2d(:))>0)
       indlayers = [indlayers i];
       props = regionprops(slice2d,'Perimeter');
       perims(i) = props(1).Perimeter*dx;
   end
end
% bottom and top
slice2d = vol_tot(:,:,indlayers(1));
props = regionprops(slice2d,'Area');
bot_area = props(1).Area*dx^2;
slice2d = vol_tot(:,:,indlayers(end));
props = regionprops(slice2d,'Area');
top_area = props(1).Area*dx^2;

surf_area = sum(perims*dzz) + bot_area + top_area;


%keyboard

%% output

data{1} = cell_length;
data{2} = cell_volume;
data{3} = single(XX);
data{4} = single(VV);
data{5} = single(tilt_tots')/pi*180;
data{6} = single(tilt_dirs')/pi*180;
data{7} = single(depths_abs_ext');
data{8} = single(depths_rel_ext');
data{9} = single(cumvols);
data{10} = single(csareas);
data{11} = surf_area;
data{12} = touch_bot;
data{13} = touch_top;
data{14} = single(anisos);
data{15} = single(ang_anisos)/pi*180;
data{16} = single(intens_alll);
data{17} = single(vs);

clear Xcell in_up intens memb qx qy qz roi0 vol vol_bulk vol_tot q_dow q_up rcell;

end

