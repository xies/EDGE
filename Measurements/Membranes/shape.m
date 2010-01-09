function [data names units] = shape(embryo, getMemb, t, z, c, dx, dz, dt, other)
% computes the tilt angles at (t, z) for Cell i given the Embryo4D,
% the membranes, and the relevant resolutions
% based on linear regresion in the vicinity of z

% the names
names{1} = 'Tilt-x';
names{2} = 'Tilt-y';
names{3} = 'Tilt-tot';
names{4} = 'Tilt-dir';
names{5} = 'All tilts-tot';
names{6} = 'All tilts-dir';
names{7} = 'Cross section area';
names{8} = 'Average cross section area';
names{9} = 'All cross section areas';
names{10} = 'Cell apical end';
names{11} = 'All fitted centroids';
names{12} = 'Depth absolute';
names{13} = 'Volume accumulate';
names{14} = 'Wedge factor'; % all, upper half, lower half
names{15} = 'Curvature';    % all, upper half, lower half
names{16} = 'Number of changes of vertex number';
names{17} = 'Average depth of vertex number changes';
names{18} = 'All numbers of vertices';
names{19} = 'All linear extensions-x';  
names{20} = 'All linear extensions-y';
names{21} = 'Anisotropy apical';  
names{22} = 'Anisotropy central'; 
names{23} = 'Anisotropy basal';  
names{24} = 'Kink depth anisotropy'; % depth absolute, relative, ...
names{25} = 'Kink slopes anisotropy'; % slope: ratio, top, bottom
names{26} = 'Anisotropy ellipse apical';  
names{27} = 'Anisotropy ellipse central'; 
names{28} = 'Anisotropy ellipse basal';  
names{29} = 'Kink depth anisotropy ellipse';
names{30} = 'Kink slopes anisotropy ellipse';
names{31} = 'Kink depth area cross section';  
names{32} = 'Kink slopes area cross section'; 
if isnan(dt)
names{33} = 'Cell length';
names{34} = 'Volume';
names{35} = 'Cell basal end';
names{36} = 'Depth relative';
names{37} = 'Cross section area central';  % central 5 micron
names{38} = 'Belly factor';
names{39} = 'Myosin depth apical';  
names{40} = 'Myosin depth relative apical'; 
names{41} = 'Myosin intensity apical'; 
names{42} = 'Myosin depth basal'; 
names{43} = 'Myosin depth relative basal'; 
names{44} = 'Myosin intensity basal'; 
end


% the units
units{1} = 'degree';
units{2} = 'degree';
units{3} = 'degree';
units{4} = 'degree';
units{5} = 'deg';
units{6} = 'deg';
units{7} = 'micron^2';
units{8} = 'micron^2';
units{9} = 'micron^2';
units{10} = 'micron';
units{11} = 'micron';
units{12} = 'micron';
units{13} = 'micron^3';
units{14} = 'micron^{-1}';
units{15} = 'micron^{-1}';
units{16} = '';
units{17} = 'micron';
units{18} = '';
units{19} = 'micron';  
units{20} = 'micron';
units{21} = ''; 
units{22} = ''; 
units{23} = ''; 
units{24} = 'micron'; 
units{25} = ''; 
units{26} = ''; 
units{27} = ''; 
units{28} = ''; 
units{29} = 'micron'; 
units{30} = ''; 
units{31} = 'micron'; 
units{32} = ''; 

data = num2cell(NaN(length(units), 1));

if isnan(dt)
units{33} = 'micron';
units{34} = 'micron^3';
units{35} = 'micron';
units{36} = '';
units{37} = 'micron^2';
units{38} = '';
units{39} = '';
units{40} = 'micron';
units{41} = '';
units{42} = '';
units{43} = 'micron';
units{44} = '';
data = num2cell(NaN(length(units), 1));
end


% if c ~= 676;% 209 %353
%  %if (c <= 358 || c>369)
%      return;
% end

%% tracked cells

% extract coordinates of centroids
centroids = Cell.centroidStack(embryo.getCellStack(c, t)) * dx;
y_values = centroids(:, 1);
x_values = centroids(:, 2);
z_values = (1:length(x_values)) * dz;
z_values = z_values(:);

% get layers that were tracked; all others are NaN in x_values 
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


%% get estimate of tilt and centroid at layer
% number of z layers chosen to calculate tilt 
dz_sel = 6;
n_sel = round(dz_sel/dz);   

% avoid layers at bottom and top 
dz_avoid = 3;
n_avoid = round(dz_avoid/dz);   

if (ntracked<n_sel+2*n_avoid);  % exit if not enough layers available
    return;
end

ind_tracked_avoid = ind_tracked(n_avoid+1:end-n_avoid);
indlb = ind_tracked_avoid(1);   % lower and 
indub = ind_tracked_avoid(end);  % lower end of bulk of cell


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
tilt_x = atand(v(1)/v(3));   % angles in x
tilt_y = atand(v(2)/v(3));   % and y
tilt_tot = atan(sqrt(v(1)^2+v(2)^2)/abs(v(3)))/pi*180;
tilt_dir = atan2(v(2),v(1))/pi*180;

% estimate centroid
Xm = mean(X,1);  % mean of centroids used for fitting
zsh = z_values(indi)-Xm(3);  % offset to current depth, in microns
xsh=zsh*v(1)/v(3);   % x and y shift in micron
ysh=zsh*v(2)/v(3);
x_est = Xm(1)+xsh;
y_est = Xm(2)+ysh;
z_est = z_values(indi);

% area of cross section
area = embryo.getCell(c, t, z).area * dx^2 * cosd(tilt_tot);


data{1} = tilt_x;
data{2} = tilt_y;
data{3} = tilt_tot;
data{4} = tilt_dir;
data{7} = area;


%% the following, do only for master layer since these properties are a 
% function of the whole cell stack, not an individual layer

if z ~= embryo.masterLayer
     return;
end

c


%% estimate tilts and centroids for all layers

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
    %keyboard
end



%% bulk of cell: calculate properties 
% (not inlcuded avoided cell at top and bottom)

avg_cs_area=0;
volume=0;
leng=0;
depths=zeros(indhl,1);
areas=zeros(indhl,1);
nverts = zeros(indhl,1);
xext = zeros(indhl,1);
yext = zeros(indhl,1);
majors = zeros(indhl,1);
minors = zeros(indhl,1);
for ind=indlb:indub   % from lower to lower part of bulk
    if tracked(ind)   % retrieve various properties
        zind = embryo.unTranslateZ(ind-1);   % index of layer
        area = embryo.getCell(c, t, zind).area * dx^2;
        verts = embryo.getCell(c, t, zind).vertexCoords * dx; 
        verts = fliplr(verts); % now the ordering is (x,y)
        [roic, R] = drawCellSmall(embryo, t, zind, c); % cell image
        xcents = Xests(ind,1:2)/dx-[R(2) R(1)];     % for elipse anisotropy
        cell_img = projimg(roic,tilt_tots(ind),tilt_dirs(ind),xcents);  
        props = regionprops(cell_img, 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
    end
    % volume area and number of vertices
    volume = volume+area*dz;         
    avg_cs_area = avg_cs_area+area*cos(tilt_tots(ind));  %cell is tilted
    leng = leng + dz/cos(tilt_tots(ind));
    depths(ind) = dz/cos(tilt_tots(ind));
    areas(ind) = area;   %areas in bulk of cell
    nverts(ind) = size(verts,1);  % number of vertices in bulk of cell
    % anisotropy based on extension in x and y
    verts3d = [verts z_values(ind)*ones(nverts(ind),1)];
    verts_trans = verts3d*0;
    for ii=1:nverts(ind)
        projv = dot(verts3d(ii,:)-Xests(ind,:),v);  % turn vertices
        verts_trans(ii,:) = verts3d(ii,:) - v'*projv; % because of tilt
    end
    [mxx,mxxi] =  max(verts_trans(:,1));  % calc max distances
    [mnx,mnxi] =  min(verts_trans(:,1));  % after projection on x, y
    xext(ind) = norm(verts_trans(mxxi,[1 3])-verts_trans(mnxi,[1 3]));
    [mxy,mxyi] =  max(verts_trans(:,2));
    [mny,mnyi] =  min(verts_trans(:,2));
    yext(ind) = norm(verts_trans(mxyi,[2 3])-verts_trans(mnyi,[2 3]));
    % anisotropy based on ellipses
    majors(ind) = props(1).MajorAxisLength*dx;
    minors(ind) = props(1).MinorAxisLength*dx;
end

% fill3(verts3d(:,2),verts3d(:,1),verts3d(:,3),0)
% hold on;
% fill3(verts_trans(:,2),verts_trans(:,1),verts_trans(:,3),0)
% hold off;
% axis equal;

%% extrapolate properties

% number of layers exrap. based on
nfit = round(8/dz*cos(mean(tilt_tots(indlb:indub))));

%areas (linear extrapoltion)
areas = extrap(indll:indhl,areas,indub,indlb,nfit);
% number of vertices 
nverts(indub+1:indhl) = nverts(indub);
nverts(indll:indlb-1) = nverts(indlb);
% x and y extension (linear extrapoltion)
xext = extrap(indll:indhl,xext,indub,indlb,nfit);
yext = extrap(indll:indhl,yext,indub,indlb,nfit);

majors = extrap(indll:indhl,majors,indub,indlb,nfit);
minors = extrap(indll:indhl,minors,indub,indlb,nfit);



%% upper end of cell: estimate length, volume, and maximal centroid

intensfac = 0.5; % factor of intensity at which extrapolation terminates
% this assumes cells have rectangular shape

% extract intensity inside cell at upper end of cell bulk
zub = embryo.unTranslateZ(indub-1);  % index of layer
intens = getMemb(t, zub);  % global memb at layer zub
siz = size(intens);
[roic, R] = drawCellSmall(embryo, t, zub, c);  % roic is cell_img
yxs = size(roic);
% shift to fitted centroid (in pixel)
R = R + round((Xests(indub,1:2)-[x_values(indub) y_values(indub)])/dx);
intensc = sum(sum(roic.*intens(R(1):R(1)+yxs(1)-1,R(2):R(2)+yxs(2)-1)));


% extract noise intensity
z_i = embryo.unTranslateZ(indhl-1);  % index of layer
intens = getMemb(t, z_i);  % global memb at layer zub
% shift R due to tilt (in pixel); R(1,2) contain y, x coordinates
Rc = R + fliplr(round((Xests(indhl,1:2)-Xests(indub,1:2))/dx));
Rc = [min([Rc(1) siz(1)-yxs(1)+1]) min([Rc(2) siz(2)-yxs(2)+1])];
Rc = [max([Rc(1) 1]) max([Rc(2) 1])];
intensc_bg = sum(sum(roic.*intens(Rc(1):Rc(1)+yxs(1)-1,Rc(2):Rc(2)+yxs(2)-1)));


% extrapolate
intensc_i = intensc;   
i=indub+1;
ints=[];
while (i<=indhl && (intensc_i-intensc_bg) >= intensfac*(intensc-intensc_bg)) 
    z_i = embryo.unTranslateZ(i-1);  % index of layer
    intens = getMemb(t, z_i);  % global memb at layer zub
    % shift R due to tilt (in pixel); R(1,2) contain y, x coordinates
    Rc = R + fliplr(round((Xests(i,1:2)-Xests(indub,1:2))/dx)); 
    Rc = [min([Rc(1) siz(1)-yxs(1)+1]) min([Rc(2) siz(2)-yxs(2)+1])];
    Rc = [max([Rc(1) 1]) max([Rc(2) 1])];
    intensc_i = sum(sum(roic.*intens(Rc(1):Rc(1)+yxs(1)-1,Rc(2):Rc(2)+yxs(2)-1)));
    depths(i) = dz/cos(tilt_tots(i));
    i=i+1;
    ints = [ints intensc_i];
end
i=i-1;

nub = (i-indub);

% cross section
area_ub = embryo.getCell(c, t, zub).area * dx^2;
avg_cs_area = avg_cs_area + (i-indub)*area_ub*cos(tilt_tots(indub));


% estimate length, volume and centroid 
length_ub = nub/cos(tilt_tots(indub))*dz;
%area_ub = embryo.getCell(c, t, zub).area * dx^2;
%volume_ub = area_ub*nub*dz;
volume_ub = sum(areas(indub+1:i))*dz;
Xub = Xests(i,:);



%% lower end of cell: estimate length, volume, and maximal centroid 

intensfac = 0.5; % factor of intensity at which extrapolation terminates
% this assumes cells have rectangular shape

% extract intensity inside cell at lower end of cell bulk
zlb = embryo.unTranslateZ(indlb-1);  % index of layer
intens = getMemb(t, zlb);  % global memb at layer zlb
siz = size(intens);
[roic, R] = drawCellSmall(embryo, t, zlb, c);  % roic is cell_img
yxs = size(roic);
% shift to fitted centroid (in pixel)
R = R + round((Xests(indlb,1:2)-[x_values(indlb) y_values(indlb)])/dx);
intensc = sum(sum(roic.*intens(R(1):R(1)+yxs(1)-1,R(2):R(2)+yxs(2)-1)));

% extract noise intensity
z_i = embryo.unTranslateZ(indll-1);  % index of layer
intens = getMemb(t, z_i);  % global memb at layer zlb
% shift R due to tilt (in pixel); R(1,2) contain y, x coordinates
Rc = R + fliplr(round((Xests(indll,1:2)-Xests(indlb,1:2))/dx));
Rc = [min([Rc(1) siz(1)-yxs(1)+1]) min([Rc(2) siz(2)-yxs(2)+1])];
Rc = [max([Rc(1) 1]) max([Rc(2) 1])];
intensc_bg = sum(sum(roic.*intens(Rc(1):Rc(1)+yxs(1)-1,Rc(2):Rc(2)+yxs(2)-1)));

if (dt~=0)
   intensc_bg = 0; 
end


% extrapolate
intensc_i = intensc;   
i=indlb-1;
ints=[];
while (i>=indll && (intensc_i-intensc_bg) >= intensfac*(intensc-intensc_bg)) 
    z_i = embryo.unTranslateZ(i-1);  % index of layer
    intens = getMemb(t, z_i);  % global memb at layer zlb
    % shift R due to tilt (in pixel); R(1,2) contain y, x coordinates
    Rc = R + fliplr(round((Xests(i,1:2)-Xests(indlb,1:2))/dx)); 
    Rc = [min([Rc(1) siz(1)-yxs(1)+1]) min([Rc(2) siz(2)-yxs(2)+1])];
    Rc = [max([Rc(1) 1]) max([Rc(2) 1])];
    intensc_i = sum(sum(roic.*intens(Rc(1):Rc(1)+yxs(1)-1,Rc(2):Rc(2)+yxs(2)-1)));
    depths(i) = dz/cos(tilt_tots(i));
    i=i-1;
    ints = [ints intensc_i];
end
i=i+1;

nlb = (indlb-i);

% cross section
area_lb = embryo.getCell(c, t, zlb).area * dx^2;
avg_cs_area = avg_cs_area + (indlb-i)*area_lb*cos(tilt_tots(indlb));

avg_cs_area = avg_cs_area/(indub-indlb+nub+nlb);  % devide by number of bulk layers

% estimate length, volume and centroid 
length_lb = nlb/cos(tilt_tots(indlb))*dz;
%area_lb = embryo.getCell(c, t, zlb).area * dx^2;
%volume_lb = area_lb*nlb*dz;
volume_lb = sum(areas(i:indlb-1))*dz;
Xlb = Xests(i,:);



%% several useful index scalars and vectors

%index of basal and apical end of cell
indbe = indlb-nlb;
indae = indub+nub;

% all indices of cell
indcell = transpose(indbe:indae);
indbulk = transpose(indlb:indub);



%% depths

% cumulative sum of distances
depths_abs = flipud(cumsum(flipud(depths)));
depths_abs = depths_abs-depths_abs(indae);
depths_abs(depths_abs<0) = 0;
depths_rel = depths_abs/max(depths_abs);

% distances extended to lowest and highest layer
depths_ext = depths;
depths_ext(indll:indbe-1)=depths(indbe);
depths_ext(indae+1:indhl)=depths(indae);
depths_abs_ext = flipud(cumsum(flipud(depths_ext)));
depths_abs_ext = depths_abs_ext-depths_abs_ext(indae);
depths_rel_ext = depths_abs_ext/max(depths_abs);

%% a few more useful indeces
indhalf = find(depths_rel<0.5);
center_ind = indhalf(1)-round(2.5/dz):indhalf(1)+round(2.5/dz);

% specific regions in cell
indbot = find(depths_rel_ext>0.7 & depths_rel_ext<0.9);
indtop = find(depths_rel_ext>0.1 & depths_rel_ext<0.3);
indmid = find(depths_rel_ext>0.4 & depths_rel_ext<0.6);


%% cumulative volume

areas_cell = areas;
areas_cell(1:indbe-1) = 0; 
areas_cell(indae+1:end) = 0;
volume_cum = flipud(cumsum(flipud(areas_cell)))*dz;
volume_cum = volume_cum-volume_cum(indae);
volume_cum(volume_cum<0) = 0;



%% cross section area
    
% cross section of central region
center_cs_area = mean(cos(tilt_tots(center_ind)) .* areas(center_ind));

areas_cs = cos(tilt_tots).*areas; % areas of all cross sections

[indkink_ac,slopetop_ac,slopebot_ac,cfun_ac,sv_ac] = ...
        fit_kink(areas_cs,indcell,indbot,indtop);
    
kinks_ac = [depths_abs_ext(indkink_ac) depths_rel_ext(indkink_ac)];
kinks_ac_slopes = [slopetop_ac/slopebot_ac slopetop_ac slopebot_ac];      


% plot for testing
%plot(indcell,areas_cs(indcell),'o',indcell,feval(cfun_ac,indcell),'-r'); 


%% wedginess and belly factor


% quantities that characterize shape of cell
wedge_fac = zeros(3,1);
wedge_fac(1) = 1+(mean(areas_cs(indbot))/mean(areas_cs(indtop))-1)/...
               abs(mean(depths_abs(indbot))-mean(depths_abs(indtop)));
wedge_fac(2) = 1+(mean(areas_cs(indmid))/mean(areas_cs(indtop))-1)/...
               abs(mean(depths_abs(indmid))-mean(depths_abs(indtop)));
wedge_fac(3) = 1+(mean(areas_cs(indbot))/mean(areas_cs(indmid))-1)/...
               abs(mean(depths_abs(indbot))-mean(depths_abs(indmid)));
belly_fac = mean(areas_cs(indmid))/mean(areas_cs([indbot;indtop]));


%% curvature of cells

% curv = zeros(indhl,1);
% for i=indbe:indae-1
%     curv(i) = norm(vs(i+1,:)-vs(i,:))/(dz/cos(tilt_tots(i)));
% end

curvat = zeros(3,1);
curvat(1) = norm(mean(vs(indbot,:),1)-mean(vs(indtop,:),1))/...
              abs(mean(depths_abs(indbot))-mean(depths_abs(indtop)));
curvat(2) = norm(mean(vs(indtop,:),1)-mean(vs(indmid,:),1))/...
              abs(mean(depths_abs(indtop))-mean(depths_abs(indmid)));
curvat(3) = norm(mean(vs(indbot,:),1)-mean(vs(indmid,:),1))/...
              abs(mean(depths_abs(indbot))-mean(depths_abs(indmid)));

          
%% vertices
nvertsshift = 2;  % >1 to avoid counting noisy fluctuating vertices
diff_verts = nverts-circshift(nverts,[nvertsshift 0]);
diff_verts([indll:nvertsshift indhl-nvertsshift:indhl]) = 0;
nverts_change = sum(abs(diff_verts))/nvertsshift;
depth_verts_change = mean(depths_abs_ext(find(diff_verts~=0)));



%% anisotropy
aniso = zeros(3,1);
aniso(1) = mean(xext(indtop))/mean(yext(indtop));
aniso(2) = mean(xext(indmid))/mean(yext(indmid));
aniso(3) = mean(xext(indbot))/mean(yext(indbot));
anisos = xext./yext;

[indkink_aniso,slopetop_aniso,slopebot_aniso,cfun_aniso,sv_aniso] = ...
          fit_kink(anisos,indcell,indbot,indtop);

kinks_aniso = [depths_abs_ext(indkink_aniso) ...
               depths_rel_ext(indkink_aniso)];
kinks_aniso_slopes = [slopetop_aniso/slopebot_aniso ...
                      slopetop_aniso slopebot_aniso];      



%% anisotropy ellipse
aniso_ellipse = zeros(3,1);
aniso_ellipse(1) = mean(majors(indtop))/mean(minors(indtop));
aniso_ellipse(2) = mean(majors(indmid))/mean(minors(indmid));
aniso_ellipse(3) = mean(majors(indbot))/mean(minors(indbot));
anisos_ellipse = majors./minors;

[indkink_aniso_ellipse,slopetop_aniso_ellipse,...
        slopebot_aniso_ellipse,cfun_aniso_ellipse, ...
               sv_aniso_ellipse] = ...
            fit_kink(anisos_ellipse,indcell,indbot,indtop);

kinks_aniso_ellipse = [depths_abs_ext(indkink_aniso_ellipse) ...
                       depths_rel_ext(indkink_aniso_ellipse)];
kinks_aniso_ellipse_slopes = ...
    [slopetop_aniso_ellipse/slopebot_aniso_ellipse ...
     slopetop_aniso_ellipse slopebot_aniso_ellipse];      
% plot for testing
%plot(indcell,anisos_ellipse(indcell),'o',indcell,feval(cfun_aniso_ellipse,indcell),'-r'); 




%% myosin at top and bottom of cell

%keyboard

if isnan(dt)
    % bottom of cell
    ints_bot=zeros(indhl,1);
    for indi=max([indbe-3 1]):min([indbe+12 indhl])
        z_i = embryo.unTranslateZ(indi-1);  % index of layer
        intens = other.YolkStalk(t, z_i);  % global YS at layer z_i
        % get closest cell outline available
        [q,ind_close] = sort(abs(indi-ind_tracked_avoid));
        ind_an = ind_tracked_avoid(ind_close(1));
        z_an = embryo.unTranslateZ(ind_an-1);
        [roic, R] = drawCellSmall(embryo, t, z_an, c);  % roic is cell_img
        yxs = size(roic);
        % shift R due to tilt (in pixel); R(1,2) contain y, x coordinates
        Rc = R + fliplr(round((Xests(indi,1:2)-Xests(ind_an,1:2))/dx));
        Rc = [min([Rc(1) siz(1)-yxs(1)+1]) min([Rc(2) siz(2)-yxs(2)+1])];
        Rc = [max([Rc(1) 1]) max([Rc(2) 1])];
        ints_bot(indi) = sum(sum(roic.*intens(Rc(1):Rc(1)+yxs(1)-1,Rc(2):Rc(2)+yxs(2)-1)));
        %Rc
    end
    ints_bot = ints_bot.*cos(tilt_tots);
    
    % top of cell
    ints_top=zeros(indhl,1);
    for indi=max([indae-12 1]):min([indae+3 indhl])
        z_i = embryo.unTranslateZ(indi-1);  % index of layer
        intens = other.YolkStalk(t, z_i);  % global YS at layer z_i
        % get closest cell outline available
        [q,ind_close] = sort(abs(indi-ind_tracked_avoid));
        ind_an = ind_tracked_avoid(ind_close(1));
        z_an = embryo.unTranslateZ(ind_an-1);
        [roic, R] = drawCellSmall(embryo, t, z_an, c);  % roic is cell_img
        yxs = size(roic);
        % shift R due to tilt (in pixel); R(1,2) contain y, x coordinates
        Rc = R + fliplr(round((Xests(indi,1:2)-Xests(ind_an,1:2))/dx));
        Rc = [min([Rc(1) siz(1)-yxs(1)+1]) min([Rc(2) siz(2)-yxs(2)+1])];
        Rc = [max([Rc(1) 1]) max([Rc(2) 1])];
        ints_top(indi) = sum(sum(roic.*intens(Rc(1):Rc(1)+yxs(1)-1,Rc(2):Rc(2)+yxs(2)-1)));
        %Rc
    end
    ints_top = ints_top.*cos(tilt_tots);

    % maximum mysoin intensities
    [mx_myo_bot,ind_myo_bot] = max(ints_bot);
    [mx_myo_top,ind_myo_top] = max(ints_top);

    myosin_depth = [depths_abs_ext(ind_myo_top) ...
                    depths_abs(1)-depths_abs_ext(ind_myo_bot)];
    myosin_depth_rel = [depths_rel_ext(ind_myo_top) ...
                        1-depths_rel_ext(ind_myo_bot)];
    myosin_intens = [mx_myo_top mx_myo_bot];
    
end

%keyboard

%% output

data{5} = tilt_tots'/pi*180;
data{6} = tilt_dirs'/pi*180;
data{8} = avg_cs_area;
data{9} = areas_cs;
data{10} = Xub;
data{11} = Xests;
data{12} = depths_abs_ext';
data{13} = volume_cum;
data{14} = wedge_fac;
data{15} = curvat;
data{16} = nverts_change;
data{17} = depth_verts_change;
data{18} = nverts;
data{19} = xext;
data{20} = yext;
data{21} = aniso(1);   % top
data{22} = aniso(2);    % mid
data{23} = aniso(3);    % bot
data{24} = kinks_aniso;
data{25} = kinks_aniso_slopes;
data{26} = aniso_ellipse(1);   % top
data{27} = aniso_ellipse(2);    % mid
data{28} = aniso_ellipse(3);    % bot
data{29} = kinks_aniso_ellipse;
data{30} = kinks_aniso_ellipse_slopes;
data{31} = kinks_ac; % area_cs_change
data{32} = kinks_ac_slopes; % 

if isnan(dt)
data{33} = leng + length_ub + length_lb;
data{34} = volume + volume_ub + volume_lb;
data{35} = Xlb;
data{36} = depths_rel_ext';
data{37} = center_cs_area;
data{38} = belly_fac;
data{39} = myosin_depth(1);
data{40} = myosin_depth_rel(1);
data{41} = myosin_intens(1);
data{42} = myosin_depth(2);
data{43} = myosin_depth_rel(2);
data{44} = myosin_intens(2);
end


%keyboard
