function mcimg = projimg(cimg,ang_tot,ang_dir,cents)
% projets cell image (cimg) on plane perpendicular to cell axis
% ang_tot is angle of cell axis with z axis
% ang_dir is angle of cell axis in xyplane 
% cents are centroids (x,y), in pixel and in coord of cimg 

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

