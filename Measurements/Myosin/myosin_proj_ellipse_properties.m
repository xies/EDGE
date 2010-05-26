function [data names units] = myosin_proj_ellipse_properties(embryo, getMyosin, t, z, c, dx, dz, dt, other)
    % computes the ellipse properties of the projected myosin at (t, z) 
    % for Cell i given the Embryo4D,
    % the membranes, and the relevant resolutions

    % the names
    names{1} = 'Ellipse x-positions';
    names{2} = 'Ellipse y-positions';
    names{3} = 'Ellipse integral';
    names{4} = 'Ellipse anisotropy';
    names{5} = 'Ellipse angle';
    names{6} = 'Ellipse relative position';
    names{7} = 'Ellipse angle of position';

    % the units
    units{1} = 'microns';
    units{2} = 'microns';
    units{3} = 'arbitrary';
    units{4} = 'dimensionless';
    units{5} = 'rad';
    units{6} = 'dimensionless';
    units{7} = 'rad';


    % for myosin projection, only the top layer needs to be processed
    if z ~= embryo.topLayer
        data = num2cell(NaN(7, 1));
        return;
    end

    myosin = getMyosin(t, z);
    
    % a few paramters
    sig_filt = 3; % sigma of kernel for smoothing myosin 
    myos_threshold = 0.8; % blobs are myosin structures above this threshold 


    [roic, R] = drawCellSmall(embryo, t, z, c);  % roic is cell_img
    % (R(1) is offset in vertical, R(2) in horizontal direction)
    yxs = size(roic);

    center = embryo.getCell(c, t, z).centroid;
    center = center(2:-1:1);

    % propty = regionprops(regions,'Centroid'); % center of cell
    % centers = [propty.Centroid];  % centers(1) is center in horizontal, centers(2) in vertical direction


% keyboard
    % extract myosin inside cell
    myosc = double(myosin(R(1):R(1)+yxs(1)-1,R(2):R(2)+yxs(2)-1));
    myosc = myosc.*roic;
    myosc = myosc/mean(mean(myosc(roic))); % normalize



    %smooth myosin & apply threshold, in order to find potential blobs 
    %and initial values of paramters for fitting:
    ext=(2^ceil(max(log2(size(myosc)))));       %smooth myosin inside cell
    [myosf,filt] = get_filtered_gauss(myosc,ext,sig_filt,0);
    [roicf,filt] = get_filtered_gauss(roic,ext,sig_filt,0);
    myosf(roic>0) = myosf(roic>0)./roicf(roic>0);
    myosf=myosf.*double(roic);

    myos_th = myosf > myos_threshold;                    %and threshold
    myosi = bwlabel(myos_th, 4);       

    %eliminate tiny blobs
    propty = regionprops(myosi,'Area');      
    areas = transpose([propty.Area]);
    qar=areas>10;
    for i=1:length(areas)
       myosi(myosi==i)=qar(i)*myosi(myosi==i); 
    end
    myos_th=myosi>0;
    myosi = bwlabel(myos_th,4);       

    %order according to total myosin intensity; put largest blob first
    nspo=max(myosi(:));                   
    sumspo=zeros(nspo,1);
    for i=1:nspo
        sumspo(i)=sum(sum(myosf(myosi==i)));
    end
    [sumspo,indspo]=sort(sumspo,'descend');
    myosi_cp=myosi;
    for i=1:nspo
        myosi(myosi_cp==indspo(i))=i;
    end



    %get initial values for fit parameters
    propty = regionprops(myosi,'Centroid');     
    cents = [propty.Centroid];
    nspots = size(cents,2)/2;
    cents = transpose(reshape(cents,[2 nspots]));
    propty = regionprops(myosi,'MajorAxisLength');
    mxax = transpose([propty.MajorAxisLength]);
    propty = regionprops(myosi,'MinorAxisLength');
    mnax = transpose([propty.MinorAxisLength]);
    propty = regionprops(myosi,'Orientation');
    orient = transpose([propty.Orientation]);
    orient = -orient/90*pi/2;


    ng=length(orient);   % number of ellipses

    xe0=1;      % set initial paramters and upper and lower bounds
    ub=3;
    lb=0;
    for i=1:ng       
        xe0=[xe0 1 cents(i,1) cents(i,2) mxax(i)/2 mnax(i)/2 orient(i)];
        ub=[ub 10 yxs(2) yxs(1) 25 25 100*2*pi];
        lb=[lb 0 1 1 0 0 -100];
    end


    if (ng > 0)   % if at least one considerable blob exists 
        opts = optimset('display', 'off');
        x = lsqcurvefit(@ellipses,xe0,roic,myosc,lb,ub, opts);   % do the fit
%         co_out=zeros(7,ng);
        co_out = cell(7,1);
        for i=1:ng              % extract paramters for all N=ng ellipses
            cf=x(2+(i-1)*6:1+i*6);
            co_out{1}(i) = (R(2)+cf(2))*dx;   % 'Myosin spot (ellipse) x position'            
            co_out{2}(i) = (R(1)+cf(3))*dx;   % 'Myosin spot (ellipse) y position'
            Fspot = ellipses([0 cf],roic);
            co_out{3}(i) = sum(Fspot(:));     % 'Myosin spot (ellipse) intensity (integral)'
            co_out{4}(i) = max([cf(4) cf(5)])/min([cf(4) cf(5)]);    % 'Myosin spot (ellipse) anisotropy (longer/shorter sigma)'
            co_out{5}(i) = cf(6);            % 'Myosin spot (ellipse) angle of ellipse'  
            [angle rel_dist] = calc_loc(roic,cents(i,:),[cf(2) cf(3)]);
            co_out{6}(i) = rel_dist;         % 'Myosin spot (ellipse) relative position in cell'
            co_out{7}(i) = angle;            % 'Myosin spot (ellipse) orientation of position'
        end
    else          
        co_out = num2cell(NaN(7,1));
    end        
%     if (ng>1)                                   % sort them according to intensity
%         [mx,indmx]=sort(co_out{3},'descend');
%         co_out=co_out(:,indmx);
%     end

    % places the data in an array
    data = co_out;

end






function F = ellipses(x,xdata)
% function fitted to myosin blobs
% sum of constants and ellipses for fitting myosin

    nx=length(x);

    a=size(xdata);
    [xk,yk] = meshgrid(1:a(2),1:a(1));

    F=xk*0+x(1);

    for i=2:6:nx
        cp=x(i);
        cx=x(i+1);
        cy=x(i+2);
        csx=x(i+3);
        csy=x(i+4);
        theta=x(i+5);
        xkt=xk-cx;
        ykt=yk-cy;
        xkturn=xkt*cos(theta)+ykt*sin(theta);
        ykturn=-xkt*sin(theta)+ykt*cos(theta);
        g=cp*exp(-(xkturn).^2/2/csx^2-(ykturn).^2/2/csy^2);
        F=F+g;
        F=F+(~xdata(round(cy),round(cx)))*10;  % penality for centering Gaussians on background
    end

    F=F.*double(xdata);


end



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