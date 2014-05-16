function [dirty,clean] = eliminate_region(regions,index2elim,borderidx)
% Eliminates the selected index2elim from the REGIONS labeled image by
% converting them to the index of the border (borderidx)

% original = regions;
for i = 1:length(index2elim)
    regions(regions == index2elim(i)) = borderidx;
end

borders = regions == borderidx;

dirty = bwmorph(borders, 'skel', Inf);
clean = bwmorph(dirty, 'shrink', Inf);
clean = bwmorph(clean, 'clean', Inf);

end
