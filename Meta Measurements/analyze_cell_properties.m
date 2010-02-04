clear all;

% parameters
dz=0.3;   % z resolution
dx=0.24;
dip=0.5;  %
xi=0:dip:100;
nip=length(xi);

xsize=150;
ysize=250;
ext=512;

times=0;


% restore measurements
dir = '../DATA_OUTPUT/';
dir_gui = '../DATA_GUI/';
%data_tags = {'b1_05z';'b1_07z';'b2_05z'};%;'b3_07z';'b3_09z';'b4_04z'};            % experiments
data_tags = {'b2_05z'};%;'b3_07z';'b3_09z';'b4_04z'};            % experiments
nd=size(data_tags,1);
% date and time measurement was taken
time_tags = 'data';  
% date and time measurement single cells were taken
time_tags_a = 'data_a';
time_tags_b = 'data_b';

% meas_tags = {'Centroid-x';'Centroid-y';'Membranes--tilt--Tilt-tot';...
%              'Nuclei--nuclei_properties--depth_relative';...
%              'Nuclei--nuclei_properties--depth';...
%              'Nuclei--nuclei_properties--length_relative';...
%              'Nuclei--nuclei_properties--length'} ;   % measured properties
meas_tags = {'Centroid-x';'Centroid-y';...
             'Membranes--shape--Tilt-x';'Membranes--shape--Tilt-y';...
             'Membranes--shape--Tilt-tot';'Membranes--shape--Tilt-dir';...
             'Membranes--vertices--Vertex-x';...
             'Membranes--vertices--Vertex-y';...
             'Membranes--vertices--# of vertices';...
             'Membranes--basic_geometry--Cell volume';...
             'Membranes--basic_geometry--Cell volume';...
             'Membranes--basic_geometry--Cell length';...
             'Membranes--basic_geometry--Cell length';...
             'Nuclei1--nuclei_properties--depth_relative';...
             'Nuclei1--nuclei_properties--length_relative';...
             'YolkStalk--yolk_stalk_properties--hole diameter';
             'Membranes--shape--Average cross section area';...
             'Membranes--shape--Cell length';...
             'Membranes--shape--Volume';
             'Membranes--shape--Cell apical end';
             'Membranes--shape--Cell basal end'};
        
ind_sel = [1 2 3 4 5 6 7 8 10 11 12 13 17 18 19 20 21];
meas_tags = meas_tags(ind_sel); 
nm=size(meas_tags,1);

%
ml=zeros(nd,1);
mlind=zeros(nd,1);
data=cell(nd,nm);
dataa=cell(nd,nm);
datab=cell(nd,nm);
embs=cell(nd,1);
for i=1:nd
    for j=1:nm
        load([dir data_tags{i} '/' time_tags '/' meas_tags{j} '.mat']);
        data{i,j}=savedata;
%         load([dir data_tags{i} '/' time_tags_a '/' meas_tags{j} '.mat']);
%         dataa{i,j}=savedata;
%         load([dir data_tags{i} '/' time_tags_b '/' meas_tags{j} '.mat']);
%         datab{i,j}=savedata;
    end
    % masterlayer
    load([dir_gui data_tags{i} '/embryo_data.mat']);
    embs{i} = embryo4d;
    ml(i) = embryo4d.masterLayer+1;
    mlind(i) = embs{i}.translateZ(ml(i))+1;
end

age = [2 3 1 1 1 2];

t=0; % time for fixed data set


%% visualization

i=1; % visualized experiment
j = find(ind_sel==18); 
meas_tags{j}  % visualized property


nc = size(data{i,1},2);  % number of cells
nz = size(data{i,1},1);  % number of layers
z_values = (1:nz)*dz;


vs = zeros(nc,1); 
for c=1:nc
    vs(c) = data{i,j}{ml(i),c}(1);
end
vs=vs';    
vss=vs;    % for color bar
vs=round(vs);
ncol=max(ceil(vs));
if (min(vs)<0)
% vs=vs-min(vs);
% vs = round(vs/max(vs)*(ncol-1)+1);

%
% color table
if (j<0)
    col=colormap(HSV(ncol));
else
    col=colormap(jet(ncol));
end


figure(2);clf;
hh=plot([0 0],[1 1],'.w'); 
vsfac = 7; % size of symbol
axis([0 ysize 0 xsize]);
%
hold on;
for c=1:nc
    if ~isnan(vs(c))
        %lay = ml(i);
        %lay = embs{i}.lowestTracked(c, t);
        lay = embs{i}.highestTracked(c, t);
        lay = embs{i}.translateZ(lay)+1;
        plot([data{i,1}{lay,c} data{i,1}{lay,c}],...
             [data{i,2}{lay,c} data{i,2}{lay,c}],...
             'o','Markersize',vsfac,...
             'MarkerFaceColor',col(vs(c),:),...
             'MarkerEdgeColor',col(vs(c),:));
    end
end
title(meas_tags{j});
colorbar;


%% Illustrating using 3D polygons top

figure(5);clf;
hold on;

%i=3;

%[min_tot,center_ind] = min(ttots);
%[ycore,xcore]=find(abs(tiltcompl)==min(min(abs(tiltcompl))));

for c=1:nc
    % select layer in each cell
    centroids = Cell.centroidStack(embs{i}.getCellStack(c, t));
    x_values = centroids(:, 2);
    indtr = find(~isnan(x_values));
    nind = length(indtr);
    lay = indtr(max([1 nind-8]));
    %lay = indtr(min([nind 10]));
    %lay = ml(i);
    xx = data{i,7}{lay,c};
    yy = data{i,8}{lay,c};
    zz = zeros(length(xx),1)+z_values(lay);
    tilttot=data{i,5}{lay,c};
    tiltdir=data{i,6}{lay,c};
    % find outermost point of polygon
    plv2d=[cosd(tiltdir) sind(tiltdir)]';
    points = [xx-xx(1) yy-yy(1)];
    dists2d = points*plv2d;
    [mx,indmx] = max(dists2d);
    xx0 = xx(indmx);
    yy0 = yy(indmx);
    % flip polygon
    plv = [plv2d*sind(tilttot);cosd(tilttot)];
    nvert = length(xx);
    coo = zeros(nvert,3);
    if ~isnan(xx)
        for v=1:nvert       
            dist = abs(dot(plv,[xx(v)-xx0 yy(v)-yy0 0]));
            coo(v,1) = xx(v)+plv(1)*dist;
            coo(v,2) = yy(v)+plv(2)*dist;
            coo(v,3) = zz(v)+plv(3)*dist;
        end
    end
    %fill3(xx,yy,zz,col(vs(c+10),:));
    if (~isnan(vs(c)))
        fill3(coo(:,2),coo(:,1),coo(:,3),col(vs(c),:));
    end
    %tilttot
    %tiltdir
end
hold off;
axis equal;
axis([0 xsize 0 ysize]);



%% Plot furrow

% vector pointing towards furrow

agecol = ['k' 'b' 'r'];

  
figure(4);clf;
hold on;
for ii=1:2:3
    nc = size(data{ii,1},2);  % number of cells
    nz = size(data{ii,i},1);  % number of layers
    z_values = (1:nz)*dz;
    cooa = zeros(3,nz);
    coob = zeros(3,nz);
    for z=1:nz
        cooa(1,z)=dataa{ii,1}{1,z};
        cooa(2,z)=dataa{ii,2}{1,z};
        cooa(3,z)=z_values(z);
        coob(1,z)=datab{ii,1}{1,z};
        coob(2,z)=datab{ii,2}{1,z};
        coob(3,z)=z_values(z);
    end
    
    inda = find(~isnan(cooa(1,:)));
    indb = find(~isnan(coob(1,:)));
    
    coa=cooa(:,inda(end));
    cob=coob(:,indb(end));
    
    fv=(coa-cob)/norm(coa-cob);
    pv=cross(fv,[0 0 1]);
    pv=pv/norm(pv);
    
    %
    
    coos=[];
    for c=1:nc
        centroids = Cell.centroidStack(embs{ii}.getCellStack(c, t));
        x_values = centroids(:, 2);
        indtr = find(~isnan(x_values));
        nind = length(indtr);
        indht = indtr(max([1 nind-5]));
        %lay = embs{ii}.highestTracked(c, t);
        %lay = ml(ii);
        %indht = embs{ii}.translateZ(lay);
        coo = [data{ii,1}{indht,c} data{ii,2}{indht,c} z_values(indht)];
        if (norm(dot(coo'-coa,pv))<4)  % distance from mid plane
            coos = [coos;coo];
        end
    end   
    coos(:,3)=coos(:,3)-coos(1,3);
    
    plot(coos(:,1),coos(:,3),[agecol(age(ii)) 'o'],'LineWidth',2,'Markersize',9)
    
end
hold off;
%axis equal;
axis([0 250 -10 30]);



%% Plot furrow, second method

% vector pointing towards furrow

agecol = ['k' 'b' 'r'];

figure(4);clf;
hold on;
for ii=2:2
    nc = size(data{ii,1},2);  % number of cells
    nz = size(data{ii,i},1);  % number of layers
    z_values = (1:nz)*dz;
    cooa = zeros(3,nz);
    coob = zeros(3,nz);
    for z=1:nz
        cooa(1,z)=dataa{ii,1}{1,z};
        cooa(2,z)=dataa{ii,2}{1,z};
        cooa(3,z)=z_values(z);
        coob(1,z)=datab{ii,1}{1,z};
        coob(2,z)=datab{ii,2}{1,z};
        coob(3,z)=z_values(z);
    end
    
    inda = find(~isnan(cooa(1,:)));
    indb = find(~isnan(coob(1,:)));
    
    coa=cooa(:,inda(end-30));
    cob=coob(:,indb(end-30));
    
    fv=(coa-cob)/norm(coa-cob);
    pv=cross(fv,[0 0 1]);
    pv=pv/norm(pv);
    
    %
    
    coos=[];
    coosl=[];
    for c=1:nc
        centroids = Cell.centroidStack(embs{ii}.getCellStack(c, t));
        x_values = centroids(:, 2);        
        indtr = find(~isnan(x_values));
        nind = length(indtr);
        dists = zeros(nind,1); 
        for zz=1:nind
           v_cent = [data{ii,1}{indtr(zz),c};data{ii,2}{indtr(zz),c};z_values(indtr(zz))];
           dists(zz) = abs(norm(cross(v_cent-coa,fv)));  
        end
        [dist inddist] = min(dists);
        if (dist<2)  % distance from mid plane
            coo = [data{ii,1}{indtr(inddist),c};...
                   data{ii,2}{indtr(inddist),c};...
                   z_values(indtr(inddist))];
            coosl = [coosl;coo'];
            coo = [data{ii,1}{indtr(nind),c};...
                   data{ii,2}{indtr(nind),c};...
                   z_values(indtr(nind))];
            coos = [coos;coo'];
            nind;
        end
    end
    coos(:,3)=coos(:,3);%-min(coos(:,3));    
    plot(coos(:,1),coos(:,3),[agecol(age(ii)) 'o'],...
        'LineWidth',2,'Markersize',9)
    coosl(:,3)=coosl(:,3);%-min(coosl(:,3));
    plot(coosl(:,1),coosl(:,3),[agecol(age(ii)) 'o'],...
        'LineWidth',2,'Markersize',2)
    
end
hold off;
%axis equal;
axis([0 250 0 50]);





%%  plot distributions of properties

agecol = ['k' 'b' 'r'];

fsz=16


figure(10);clf;


jj = find(ind_sel==10); %volume

hold on;
for ii=1:nd
    nc = size(data{ii,1},2);  % number of cells
    vs = zeros(nc,1); 
    for c=1:nc
        vs(c) = data{ii,jj}{ml(ii),c};  % volume # 10
    end
    plot([sort(vs)' 1200],[(1:nc)/nc 1],agecol(age(ii)),'LineWidth',2)
end
hold off;
axis([0 1100 0 1])
xlabel('Volume  [\mum^3]','fontsize',fsz)
ylabel('Cumulative frequency','fontsize',fsz)
set(gca,'fontsize',fsz);
legend(char(data_tags),'Location','NorthWest');

%%  
figure(11);clf;

jj = find(ind_sel==11); %length

hold on;
for ii=1:nd
    nc = size(data{ii,1},2);  % number of cells
    vs = zeros(nc,1); 
    for c=1:nc
        vs(c) = data{ii,jj}{ml(ii)+1,c};  % volume # 10
    end
    plot([sort(vs)' 50],[(1:nc)/nc 1],agecol(age(ii)),'LineWidth',2)
end
hold off;

axis([0 50 0 1])
xlabel('Length  [\mum]','fontsize',fsz)
ylabel('Cumulative frequency','fontsize',fsz)
set(gca,'fontsize',fsz);
legend(char(data_tags),'Location','NorthWest');



%%  
figure(12);clf;

jj = find(ind_sel==12); %nuclei depth

hold on;
for ii=1:nd
    nc = size(data{ii,1},2);  % number of cells
    vs = zeros(nc,1); 
    for c=1:nc
        vs(c) = data{ii,jj}{ml(ii)+1,c};  % volume # 10
    end
    plot([sort(vs)' 50],[(1:nc)/nc 1],agecol(age(ii)),'LineWidth',2)
end
hold off;

axis([0 1 0 1])
xlabel('Nuclei depth (relative to cell length)','fontsize',fsz)
ylabel('Cumulative frequency','fontsize',fsz)
set(gca,'fontsize',fsz);
legend(char(data_tags),'Location','NorthWest');



%%  
figure(13);clf;

jj = find(ind_sel==13); %nuclei length

hold on;
for ii=1:nd
    nc = size(data{ii,1},2);  % number of cells
    vs = zeros(nc,1); 
    for c=1:nc
        vs(c) = data{ii,jj}{ml(ii)+1,c};  % volume # 10
    end
    plot([sort(vs)' 50],[(1:nc)/nc 1],agecol(age(ii)),'LineWidth',2)
end
hold off;

axis([0 1 0 1])
xlabel('Nuclei length (relative to cell length)','fontsize',fsz)
ylabel('Cumulative frequency','fontsize',fsz)
set(gca,'fontsize',fsz);
legend(char(data_tags),'Location','NorthWest');



%%  
figure(14);clf;

jj = find(ind_sel==14); %nuclei depth

hold on;
for ii=1:6
    nc = size(data{ii,1},2);  % number of cells
    vs = zeros(nc,1); 
    for c=1:nc
        vs(c) = data{ii,jj}{ml(ii)+1,c};  % volume # 10
    end
    plot([sort(vs)' 50],[(1:nc)/nc 1],agecol(age(ii)),'LineWidth',2)
end
hold off;

axis([0 3 0 1])
xlabel('Yolk Stalk Radius [microns]','fontsize',fsz)
ylabel('Cumulative frequency','fontsize',fsz)
set(gca,'fontsize',fsz);
legend(char(data_tags),'Location','NorthWest');










%% plot contour angles



%% visualization

i=1; % visualized experiment
%j=11; % visualized property
j = find(ind_sel==18); %volume
meas_tags{j}

% if (j==5 || j==6)
%     ml(1)=50;
% end

nc = size(data{i,1},2);  % number of cells
nz = size(data{i,1},1);  % number of layers
z_values = (1:nz)*dz;



% smoothed tilt for tilt contour lines
smfac=25;
tiltx=zeros(xsize,ysize);
tilty=zeros(xsize,ysize);
tilttot=zeros(xsize,ysize);
tiltdir=zeros(xsize,ysize);
xxs=zeros(nc,1);
yys=zeros(nc,1);
ttots = zeros(nc,1);

for c=1:nc
    centroids = Cell.centroidStack(embs{i}.getCellStack(c, t));
    x_values = centroids(:, 2);
    indtr = find(~isnan(x_values));
    nind = length(indtr);
    %lay = embs{i}.highestTracked(c, t);            
    lay = mlind(i);
    %lay = embs{i}.translateZ(lay)+1;
    %lay = indtr(max([1 nind-10]));
    xx = data{i,1}{lay,c};
    yy = data{i,2}{lay,c};
    tiltx(round(yy),round(xx))=data{i,3}{lay,c};
    tilty(round(yy),round(xx))=data{i,4}{lay,c};
    tilttot(round(yy),round(xx))=data{i,5}{lay,c};
    tiltdir(round(yy),round(xx))=data{i,6}{lay,c};
    xxs(c)=xx;
    yys(c)=yy;
    ttots(c) = data{i,5}{lay,c};
end
tilty(isnan(tilty))=0;
tiltx(isnan(tiltx))=0;

indn=tiltx*0;
indn(tiltx~=0)=1;
[tiltxf,filt] = get_filtered_gauss(tiltx,ext,smfac,0);
[tiltyf,filt] = get_filtered_gauss(tilty,ext,smfac,0);
[tilttotf,filt] = get_filtered_gauss(tilttot,ext,smfac,0);
[tiltdirf,filt] = get_filtered_gauss(tiltdir,ext,smfac,0);
[tiltcompl,filt] = get_filtered_gauss_compl(tilttot.*...
                    exp(complex(0,1)*tiltdir/180*pi),ext,smfac,0);
[indnf,filt] = get_filtered_gauss(indn,ext,smfac,0);
tiltxf(indnf>0)=tiltxf(indnf>0)./indnf(indnf>0);
tiltyf(indnf>0)=tiltyf(indnf>0)./indnf(indnf>0);
tilttotf(indnf>0)=tilttotf(indnf>0)./indnf(indnf>0);
tiltdirf(indnf>0)=tiltdirf(indnf>0)./indnf(indnf>0);
tiltcompl(indnf>0)=tiltcompl(indnf>0)./indnf(indnf>0);
[roi,filt] = get_filtered_gauss(indn,ext,5,0);
roi=roi>0.1*max(roi(:));

% visualized property
vs = zeros(nc,1); 
for c=1:nc
    vs(c) = data{i,j}{ml(i)+1,c}(1);
end
vs=vs';    
vss=vs;    % for color bar
vs=round(vs);
ncol=max(ceil(vs));
% vs=vs-min(vs);
% vs = round(vs/max(vs)*(ncol-1)+1);

%
% color table
if (j<0)
    col=colormap(HSV(ncol));
else
    col=colormap(jet(ncol));
end







figure(2);clf;
hh=plot([0 0],[1 1],'.w'); 
vsfac = 7; % size of symbol
axis([0 ysize 0 xsize]);
%
hold on;
for c=1:nc
    if ~isnan(vs(c))
        %lay = ml(i);
        %lay = embs{i}.lowestTracked(c, t);
        lay = embs{i}.highestTracked(c, t);
        lay = embs{i}.translateZ(lay)+1;
        plot([data{i,1}{lay,c} data{i,1}{lay,c}],...
             [data{i,2}{lay,c} data{i,2}{lay,c}],...
             'o','Markersize',vsfac,...
             'MarkerFaceColor',col(vs(c),:),...
             'MarkerEdgeColor',col(vs(c),:));
    end
end
title(meas_tags{j});
% countour plot for tilts
[Yc,Xc] = meshgrid(min(round(xxs)):max(round(xxs)),...
                   min(round(yys)):max(round(yys)));
               
[C,h] = contour(Yc,Xc,tiltxf(min(round(yys)):max(round(yys)),...
                   min(round(xxs)):max(round(xxs))),'k','LineWidth',2);
clabel(C,h,'FontSize',15,'Rotation',0);
[C,h] = contour(Yc,Xc,tiltyf(min(round(yys)):max(round(yys)),...
                  min(round(xxs)):max(round(xxs))),'k','LineWidth',2);
clabel(C,h,'FontSize',15,'Rotation',0);


hold off;

%

figure(3);clf;
hh=plot([0 0],[1 1],'.w'); 
vsfac = 7; % size of symbol
axis([0 ysize 0 xsize]);
%
hold on;
for c=1:nc
    if ~isnan(vs(c))
        %lay = ml(i);
        %lay = embs{i}.lowestTracked(c, t);
        lay = embs{i}.highestTracked(c, t);
        lay = embs{i}.translateZ(lay)+1;
        plot([data{i,1}{lay,c} data{i,1}{lay,c}],...
             [data{i,2}{lay,c} data{i,2}{lay,c}],...
             'o','Markersize',vsfac,...
             'MarkerFaceColor',col(vs(c),:),...
             'MarkerEdgeColor',col(vs(c),:));
    end
end
title(meas_tags{j});
% countour plot for tilts
[Yc,Xc] = meshgrid(min(round(xxs)):max(round(xxs)),...
                   min(round(yys)):max(round(yys)));
               
[C,h] = contour(Yc,Xc,abs(tiltcompl(min(round(yys)):max(round(yys)),...
                   min(round(xxs)):max(round(xxs)))),'k','LineWidth',2);
clabel(C,h,'FontSize',15,'Rotation',0);
%  [C,h] = contour(Yc,Xc,180/pi*angle(tiltcompl(min(round(yys)):max(round(yys)),...
%                    min(round(xxs)):max(round(xxs)))),'k','LineWidth',2);
%  clabel(C,h,'FontSize',15,'Rotation',0);
















%%

figure(11);clf;
hold on;
for ii=nd:-1:1
    leng=length(lens{i}(sels{i}==1));
    plot([sort(lens{i}(sels{i}==1))' 50],[(1:leng)/leng 1],'LineWidth',2,'Color',colo(i,:))
end
hold off;
axis([0 50 0 1])
xlabel('Length  [\mum]','fontsize',fsz)
ylabel('Cumulative frequency','fontsize',fsz)
set(gca,'fontsize',fsz);
legend(num2str(nsel),'Location','SouthEast')
%legend('06 blast','15 blast','28 blast','23 early gs','08 early gs','21 mid gs','19 mid gs','31 mid gs','13 mid gs','04 late gs','Location','SouthEast');



%%  
figure(11);clf;

jj = find(ind_sel==11); %length

hold on;
for ii=1:nd
    nc = size(data{ii,1},2);  % number of cells
    vs = zeros(nc,1); 
    for c=1:nc
        vs(c) = data{ii,jj}{ml(ii)+1,c};  % volume # 10
    end
    plot([sort(vs)' 50],[(1:nc)/nc 1],agecol(age(ii)),'LineWidth',2)
end
hold off;

axis([0 50 0 1])
xlabel('Length  [\mum]','fontsize',fsz)
ylabel('Cumulative frequency','fontsize',fsz)
set(gca,'fontsize',fsz);
legend(char(data_tags),'Location','NorthWest');





%%

% furrow region
yrl=[60 60 60   60 55   85 75 75 75   80];
yru=[110 100 110   100 95   105 90 90 90  100];
xrl=[20 20 40   10 20   0 0 0 0  40];
xru=[220 220 220   220 230   250 250 250 250   220];


% master layer
lay=[32 46 38   36 33   29 42 44 55   30];

age=[0 0 0  2 2  4 4 4 4  6];
age=age/max(age);

colo=[0.3 0.3 0.3 ; 0.0 0.0 0.0 ; 0.6 0.6 0.6 ; ...
      0.0 0.7 1.0 ; 0.0 0.0 1.0 ; ...
      0.8 1.0 0.0 ; 0.6 1.0 0.2 ; 1.0 1.0 0.0 ; 0.4 1.0 0.4 ; ...
      1.0 0.0 0.0];

xfu=[62 122 66 47 118 123 120];
yfu=[43 93 65 44 86 83 103];
xff=[98 176 86 95 175 173 140];
yff=[43 96 68 48 81 83 103];


%%


%%
fsz=16;

i1=1;
i2=10;

mmar=cell(nd,2);
mmwx=cell(nd,2);
mmwy=cell(nd,2);
mmxy=cell(nd,2);
mmnv=cell(nd,2);
ths=[20 20 20   20 20    20 20 20 20 20 ];

for i=i1:i2
    ipar=ipars{i};
    ipwx=ipwxs{i};
    ipwy=ipwys{i};
    ipnv=ipnvs{i};
    sel=sels{i};
    qmm=zeros(nip,2);
    qwx=zeros(nip,2);
    qwy=zeros(nip,2);
    qxy=zeros(nip,2);
    qnv=zeros(nip,2);
    for j=1:nip
       rr=find(sel==1 & ~isnan(ipar(:,j)));
       if (length(rr) > ths(i))
           qmm(j,1)=mean(ipar(rr,j)) ;
           qmm(j,2)=std(ipar(rr,j)); 
           qwx(j,1)=mean(ipwx(rr,j)) ;
           qwx(j,2)=std(ipwx(rr,j)); 
           qwy(j,1)=mean(ipwy(rr,j)) ;
           qwy(j,2)=std(ipwy(rr,j)); 
           qxy(j,1)=mean(ipwx(rr,j)./ipwy(rr,j)) ;
           qxy(j,2)=std(ipwx(rr,j)./ipwy(rr,j)); 
           qnv(j,1)=mean(ipnv(rr,j)) ;
           qnv(j,2)=std(ipnv(rr,j)); 
       else
           qmm(j,1)=NaN ;
           qmm(j,2)=NaN;            
           qwx(j,1)=NaN;
           qwx(j,2)=NaN;            
           qwy(j,1)=NaN;
           qwy(j,2)=NaN;            
           qxy(j,1)=NaN;
           qxy(j,2)=NaN;            
           qnv(j,1)=NaN;
           qnv(j,2)=NaN;            
       end
    end
    mmar{i}=qmm;
    mmwx{i}=qwx;
    mmwy{i}=qwy;
    mmxy{i}=qxy;
    mmnv{i}=qnv;
end
%
figure(12);clf;
hold on;
q=ipar.*repmat(sel,[1 nip]);
plot(q')
hold off;
axis([0 60 0 50])

figure(13);clf;
hold on;
for i=i1:i2    
    plot(xi,mmar{i}(:,1),'LineWidth',3,'Color',colo(i,:));    
    plot(xi,mmar{i}(:,2),'--','LineWidth',1.0,'Color',colo(i,:));    
end
hold off;
axis([-2 40 0 60])

xlabel('Distance from apical end  [\mum]','fontsize',fsz)
ylabel('Area cross-section  [\mum^2]','fontsize',fsz)
%legend(','Location','SouthEast')


figure(14);clf;
hold on;
for i=i1:i2    
    plot(xi,mmwx{i}(:,1),'LineWidth',3,'Color',colo(i,:));    
    plot(xi,mmwx{i}(:,2),'--','LineWidth',1.0,'Color',colo(i,:));    
end
hold off;
axis([-2 40 0 10])

xlabel('Distance from apical end  [\mum]','fontsize',fsz)
ylabel('Linear extension (AP)  [\mum]','fontsize',fsz)


figure(15);clf;
hold on;
for i=i1:i2    
    plot(xi,mmwy{i}(:,1),'LineWidth',3,'Color',colo(i,:));    
    plot(xi,mmwy{i}(:,2),'--','LineWidth',1.0,'Color',colo(i,:));    
end
hold off;
axis([-2 40 0 10])

xlabel('Distance from apical end  [\mum]','fontsize',fsz)
ylabel('Linear extension (VD)  [\mum]','fontsize',fsz)


figure(16);clf;
hold on;
for i=i1:i2    
    plot(xi,mmxy{i}(:,1),'LineWidth',3,'Color',colo(i,:));    
    plot(xi,mmxy{i}(:,2),'--','LineWidth',1.0,'Color',colo(i,:));    
end
hold off;
axis([-2 40 0 3])

xlabel('Distance from apical end  [\mum]','fontsize',fsz)
ylabel('Anisotropy (AP / VD)','fontsize',fsz)


figure(17);clf;
hold on;
for i=i1:i2    
    plot(xi,mmnv{i}(:,1),'LineWidth',3,'Color',colo(i,:));    
    plot(xi,mmnv{i}(:,2),'--','LineWidth',1.0,'Color',colo(i,:));    
end
hold off;
axis([-2 40 0 8])

xlabel('Distance from apical end  [\mum]','fontsize',fsz)
ylabel('Number vertices  [\mum]','fontsize',fsz)

%%

nbins=5;
xa=zeros(nbins*2,nd);
ya=zeros(nbins*2,nd);
for j=1:nd
   
   qx=xxs{j}(sels{j}==1)-xfu(j);
   qy=yys{j}(sels{j}==1)-yfu(j);
   qa=vols{j}(sels{j}==1); 
   [qxa,xb] = bin_vol(qx,qa,nbins); 
   xa(:,j)=qxa;
   [qya,yb] = bin_vol(qy,qa,nbins); 
   ya(:,j)=qya;
end

figure(8);clf;
subplot(1,2,1);
plot(xb,xa(:,1));
hold on;
for i=1:nd
    plot(xb,xa(:,i),'LineWidth',3,'Color',[0 age(i) age(i)]);
end
axis([-100 100 400 700])
title('xdep')
hold off

subplot(1,2,2);
plot(yb,ya(:,1));
hold on;
for i=1:nd
    plot(yb,ya(:,i),'LineWidth',3,'Color',[0 age(i) age(i)]);
end
axis([-40 40 400 700])
title('ydep')
hold off
