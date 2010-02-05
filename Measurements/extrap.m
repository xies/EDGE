function yfit = extrap(ind,y,indub,indlb,nfit)
% extraplates properties from bulk region to parts 
% below and above bulk

indhl = ind(end);
indll = ind(1);

yfit = y;

% upper part
indfitu = indub-nfit+1:indub;
pfit = polyfit(indfitu',y(indfitu),1);
yfit(indub+1:indhl) = polyval(pfit,indub+1:indhl);
% lower part
indfitl = indlb:indlb+nfit-1;
pfit = polyfit(indfitl',y(indfitl),1);
yfit(indll:indlb-1) = polyval(pfit,indll:indlb-1);


end

