% Change the extension of the filename string "in". If "type" is 'none',
% then changes to no extension. If "type" is 'mat', changes to .mat.
% Returns the new filename.

function out = chgext(in, type)

% remove the extension from "in"
dot_in_string = strfind(in, '.');
dot_in_string = dot_in_string(end); % just in case 'Matlab' appears upstream in the path
noext = in(1:dot_in_string-1);

if strcmp(type, 'none')
    out = noext;
elseif strcmp(type, 'mat')
    out = strcat(noext, '.mat');
else % should not happen
    out = in;
end