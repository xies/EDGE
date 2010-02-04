clear all;

% parameters
dz=0.3;   % z resolution
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
data_tags = {'b2_05z'};            % experiments
nd=size(data_tags,1);
time_tags = {'12-6-2009, 4;19 PM'};  % date and time measurement was taken
% meas_tags = {'Centroid-x';'Centroid-y';'Membranes--tilt--Tilt-tot';...
%              'Nuclei--nuclei_properties--depth_relative';...
%              'Nuclei--nuclei_properties--depth';...
%              'Nuclei--nuclei_properties--length_relative';...
%              'Nuclei--nuclei_properties--length'} ;   % measured properties
meas_tags = {'Centroid-x';'Centroid-y';...
             'Membranes--tilt--Tilt-x';'Membranes--tilt--Tilt-y';...
             'YolkStalk--yolk_stalk_properties--x-positions';...
             'YolkStalk--yolk_stalk_properties--y-positions';...
             'YolkStalk--yolk_stalk_properties--angle';...
             'YolkStalk--yolk_stalk_properties--anisotropy';...
             'YolkStalk--yolk_stalk_properties--hole diameter';...
             'YolkStalk--yolk_stalk_properties--intensity';...
             'YolkStalk--yolk_stalk_properties--volume outside cell'};
nm=size(meas_tags,1);

%
ml=zeros(nd);
data=cell(nd,nm);
embs=cell(nd,1);
for i=1:nd
    for j=1:nm
        load([dir data_tags{i} '/' time_tags{i} '/' meas_tags{j} '.mat']);
        data{i,j}=savedata;
    end
    % masterlayer
    load([dir_gui data_tags{i} '/embryo_data.mat']);
    ml(i) = embryo4d.masterLayer;
    embs{i} = embryo4d;
end



%% visualization

i=1; % visualized experiment
j=6; % visualized property

% if (j==5 || j==6)
%     ml(1)=50;
% end

t=0; % time for fixed data set

nc = size(data{1,1},2);  % number of cells

% smoothed tilt for tilt contour lines
smfac=25;
tiltx=zeros(xsize,ysize);
tilty=zeros(xsize,ysize);
xxs=zeros(nc,1);
yys=zeros(nc,1);
for c=1:nc
    %lay = embs{i}.highestTracked(c, t);            
    lay = ml(i);
    lay = embs{i}.translateZ(lay)+1;
    xx = data{i,1}{lay,c};
    yy = data{i,2}{lay,c};
    tiltx(round(yy),round(xx))=data{i,3}{lay,c};
    tilty(round(yy),round(xx))=data{i,4}{lay,c};
    xxs(c)=xx;
    yys(c)=yy;
end
tilty(isnan(tilty))=0;
tiltx(isnan(tiltx))=0;

indn=tiltx*0;
indn(tiltx~=0)=1;
[tiltxf,filt] = get_filtered_gauss(tiltx,ext,smfac,0);
[tiltyf,filt] = get_filtered_gauss(tilty,ext,smfac,0);
[indnf,filt] = get_filtered_gauss(indn,ext,smfac,0);
tiltxf(indnf>0)=tiltxf(indnf>0)./indnf(indnf>0);
tiltyf(indnf>0)=tiltyf(indnf>0)./indnf(indnf>0);
[roi,filt] = get_filtered_gauss(indn,ext,5,0);
roi=roi>0.1*max(roi(:));

%
ncol=64;

if (j==9)
    col=colormap(HSV(ncol));
else
    col=colormap(jet(ncol));
end

% visualized property
vs = zeros(nc,1); 
for c=1:nc
    vs(c) = data{i,j}{ml(i)+1,c};
end
vs=vs';    
vs=vs-min(vs);
vs = round(vs/max(vs)*(ncol-1)+1);

figure(1);clf;
plot([0 0],[1 1],'.w'); 
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
             'MarkerEdgeColor',col(vs(c),:))
    end
end

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




%%




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


%calculate properties

nars=cell(nd,1);
vols=cell(nd,1);
sels=cell(nd,1);
xxs=cell(nd,1);
yys=cell(nd,1);
ars=cell(nd,1);
aars=cell(nd,1);
dzs=cell(nd,1);
lens=cell(nd,1);
ipars=cell(nd,1);
ipwxs=cell(nd,1);
ipwys=cell(nd,1);
ipnvs=cell(nd,1);
nsel=zeros(nd,1);
for j=1:nd
    data=dataa{j};
    dxc=datax{j};
    dyc=datay{j};
    tiltx=datatiltx{j};
    tilty=datatilty{j};
    widx=datawidx{j};
    widy=datawidy{j};
    nvert=datanvert{j};
    n=length(data);
    nn=length(data{1}.single);
    
 %   areas=zeros(n,nn);
    nar=zeros(n,1);
    vol=zeros(n,1);
    sel=zeros(n,1);
    xx=zeros(n,1);
    yy=zeros(n,1);
    ddz=zeros(n,1);    
    len=zeros(n,1);    
    aar=zeros(n,nn);
    ipar=zeros(n,nip);    
    ipwx=zeros(n,nip);    
    ipwy=zeros(n,nip);    
    ipnv=zeros(n,nip);    
    for i=1:n       
        ar=data{i}.single';
        wx=widx{i}.single';
        wy=widy{i}.single';
        nv=nvert{i}.single';
        tiltfac=cosd(tiltx{i}.single(1))*cosd(tilty{i}.single(1));
        aar(i,:)=ar*tiltfac;
        ddz(i)=dz/tiltfac;
        %if (j==3) 
        %    ar=ar*2;
        %end
        mn=min(find(~isnan(ar)));
        mx=max(find(~isnan(ar)));
        inda=find(~isnan(ar));
        nar(i)=mx-mn+1;
        len(i)=nar(i)*dz/tiltfac;
        %interpolation 
        if (nar(i)>5)
            x=(max(inda)-flipud(inda))*dz/tiltfac;
            % area
            Y=flipud(ar(inda))*tiltfac;  
            yi = interp1(x(~isnan(x)),Y(~isnan(Y)),xi);
            ipar(i,:)=yi;
            % width-x
            Y=flipud(wx(inda))*cosd(tiltx{i}.single(1));
            yi = interp1(x(~isnan(x)),Y(~isnan(Y)),xi);
            ipwx(i,:)=yi;
            % width-y
            Y=flipud(wy(inda))*cosd(tilty{i}.single(1));
            yi = interp1(x(~isnan(x)),Y(~isnan(Y)),xi);
            ipwy(i,:)=yi;
            % numb vertices
            Y=flipud(nv(inda));
            yi = interp1(x(~isnan(x)),Y(~isnan(Y)),xi);
            ipnv(i,:)=yi;            
        end
        %volume
        volo=0;
        for l=mn:mx
            if (~isnan(ar(l)))
                volo=volo+ar(l);                
            else
                q=inda(find(inda<l));
                indpre=q(end);
                q=inda(find(inda>l));
                indpost=q(1);
                volo=volo+(ar(indpre)+ar(indpost))/2;
               if isnan(volo) 
                   keyboard
               end
            end
        end
        vol(i)=volo*dz;   
%        keyboard
        xx(i)=dxc{i}.single(lay(j));
        yy(i)=dyc{i}.single(lay(j));
        sel(i) = (xx(i)>xrl(j) & xx(i)<xru(j) & ...
                    yy(i)>yrl(j) & yy(i)<yru(j) & ...
                     nar(i)>20 & vol(i)>200);   
        ars{j}(i)=ar(lay(j));
    end
    
    nars{j} = nar;
    vols{j} = vol;
    sels{j} = sel;
    xxs{j} = xx;
    yys{j} = yy;
    aars{j} = aar;
    dzs{j} = ddz;
    lens{j} = len;
    ipars{j} = ipar;
    ipwxs{j} = ipwx;
    ipwys{j} = ipwy;
    ipnvs{j} = ipnv;
    nsel(j) = sum(sels{j});
end

nsel
%xxs{4}=-xxs{4}+250


%%
% plot volumes, number of slices, area at masterlayer; 
% cells selected and in predefined bounds: red

j=5;
vsfac=10;
figure(7);clf;
subplot(1,3,1);
plot(xxs{j},yys{j},'.w')
hold on;
vs=vols{j}/max(vols{j});   % volume
for i=1:length(xxs{j})
    plot([xxs{j}(i) xxs{j}(i)],[yys{j}(i) yys{j}(i)],'o',...
        'Markersize',vs(i)*vsfac,'MarkerFaceColor',[vs(i) 0 0],'MarkerEdgeColor',[vs(i) 0 0])
end
for i=1:length(xxs{j})
    if (sels{j}(i)==0)
        plot([xxs{j}(i) xxs{j}(i)],[yys{j}(i) yys{j}(i)],'o',...
            'Markersize',vs(i)*vsfac,'MarkerFaceColor',[0 vs(i) 0],'MarkerEdgeColor',[0 vs(i) 0])
    end
end
hold off;
title('volumes')

subplot(1,3,2);
plot(xxs{j},yys{j},'.w')
hold on;
vs=nars{j}/max(nars{j});  % number of slices
for i=1:length(xxs{j})
    plot([xxs{j}(i) xxs{j}(i)],[yys{j}(i) yys{j}(i)],'o',...
        'Markersize',vs(i)*vsfac,'MarkerFaceColor',[vs(i) 0 0],'MarkerEdgeColor',[vs(i) 0 0])
end
for i=1:length(xxs{j})
    if (sels{j}(i)==0)
        plot([xxs{j}(i) xxs{j}(i)],[yys{j}(i) yys{j}(i)],'o',...
            'Markersize',vs(i)*vsfac,'MarkerFaceColor',[0 vs(i) 0],'MarkerEdgeColor',[0 vs(i) 0])
    end
end
hold off;
title('number of slices')

subplot(1,3,3);
plot(xxs{j},yys{j},'.w')
hold on;
vs=ars{j}/max(ars{j});  % area at masterlayer
for i=1:length(xxs{j})
    if (~isnan(ars{j}(i)))
    plot([xxs{j}(i) xxs{j}(i)],[yys{j}(i) yys{j}(i)],'o',...
        'Markersize',vs(i)*vsfac,'MarkerFaceColor',[vs(i) 0 0],'MarkerEdgeColor',[vs(i) 0 0])
    end
end
hold off;
title('area at master depth')



%%

i1=8;
i2=8;

figure(1);clf;
%plot(vols{2}.*sels{2},'*g')
hold on;    
for i=i1:i2
    plot(vols{i}.*sels{i},'.','Color',colo(i,:))
end
hold off;
title('volumes')

figure(2);clf;
hold on;    
for i=i1:i2
    plot(nars{i}.*sels{i},'.','Color',colo(i,:))
end
hold off;
title('number of slices')

%
figure(3);clf;
hold on;    
for i=i1:i2
    plot(nars{i}.*sels{i},vols{i}.*sels{i},'.','Color',colo(i,:))
end
hold off;
%title('volumes')


figure(4);clf;
hold on;    
for i=i1:i2
    plot(ars{i}'.*sels{i},'.','Color',colo(i,:))
end
hold off;

title('area at master depth')


figure(5);clf;
hold on;   
for i=i1:i2
    plot(xxs{i}.*sels{i},vols{i}.*sels{i},'.','Color',colo(i,:))
end
hold off;

title('x')

figure(6);clf;
hold on;    
for i=i1:i2
    plot(yys{i}.*sels{i},vols{i}.*sels{i},'.','Color',colo(i,:))
end
hold off;

title('y')

figure(8);clf;
hold on;   
for i=i1:i2
    plot(xxs{i}.*sels{i},vols{i}.*sels{i},'.','Color',colo(i,:))
end
hold off;

title('x')

figure(9);clf;
hold on;   
for i=i1:i2
    plot(lens{i}.*sels{i},'.','Color',colo(i,:))
end
hold off;

title('lengths')



%%

figure(10);clf;
fsz=16

hold on;
for i=1:nd
    leng=length(vols{i}(sels{i}==1));
    plot([sort(vols{i}(sels{i}==1))' 1200],[(1:leng)/leng 1],'LineWidth',2,'Color',colo(i,:))
end
hold off;
axis([0 1100 0 1])
xlabel('Volume  [\mum^3]','fontsize',fsz)
ylabel('Cumulative frequency','fontsize',fsz)
set(gca,'fontsize',fsz);
legend('06 blastoderm','15 blastoderm','28 blastoderm','23 early gastrula','08 early gastrula',...
        '21 mid gastrula','19 mid gastrula','31 mid gastrula','13 mid gastrula','04 late gastrula',...
            'Location','NorthWest');


figure(11);clf;
hold on;
for i=nd:-1:1
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
