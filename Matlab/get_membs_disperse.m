function [mask,status_mse,status_skel] = get_membs_disperse(filename, cutoff, X,Y)


cutoff = 10;

[file_dir,name,ext] = fileparts( filename );
im = imread(filename);
im = double(im); im = im / max(max(im));
imwrite(im,filename);

Nsmooth = 3;
% Calculate Morse Smale complex, cut at cutoff, and convert the filament structure:

mse_cmds = strjoin({ ...
    'mse' , filename, ...
    '-outDir', file_dir, ...
    '-cut', num2str(cutoff), ...
    ' -dumpArcs CUD'} , ' ');

[status_mse,result] = system(mse_cmds);

%
skl_filename = [filename '_c' num2str(cutoff) '.CUD.NDskl'];

skl_cmds = strjoin({ ...
    'skelconv', skl_filename, ...
    '-outDir', file_dir, ...
    '-smooth', num2str(Nsmooth), ...
    '-to', 'NDskl_ascii'}, ' ');

[status_skel,result] = system(skl_cmds);

output_filename = [filename ...
    '_c' num2str(cutoff) ...
    '.CUD.NDskl.S003.a.NDskl'];

MSC = parseNDskl(output_filename);

mask = MSC2mask(MSC, X,Y);
figure,imshow(imoverlay(im,mask,[1 0 0]))
% keyboard
% save([file_dir, '/MSC.c_' num2str(cutoff) '.CUD.NDskl.S003.a.NDskl.mat'],'MSC');

end