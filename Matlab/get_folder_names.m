function out = get_folder_names(src)
% gets the names of the datasets by checking the names of the folders
% in the directory src. returns them in a cell array, which is then used
% to create the dropdown menu

files = dir(src);  % like pwd, but much better for this application

goodfiles = zeros(length(files), 1);
% get names of the directories
for i = 1:length(files)    
    if files(i).isdir && ~strcmp(files(i).name , '.') && ~strcmp(files(i).name, '..')
        % if it's a directory but not the '.' directory
        goodfiles(i) = 1;
    end
end
goodfiles = find(goodfiles);

out = cell(length(goodfiles), 1);
for i = 1:length(out)
    out{i} = files(goodfiles(i)).name;
end


        