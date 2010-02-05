function [indkink,slopetop,slopebot,cfun,start_values] = fit_kink(y,indcell,indbot,indtop)
% fits a piecewise quasi linear model to a 1D cell property. 

% used fitting function
funct = 'c1+b1*((x-a1)-sqrt(1+(x-a1).^2))+b2*((x-a1)+sqrt(1+(x-a1).^2))';
fit_model = fittype(funct);   % create fit model
coeffnames(fit_model);   % returns order of paramters 

% priming paramters
pfit = polyfit(indbot,y(indbot),1);
b1i = pfit(1);                 % slope at lower part of cell
pfit = polyfit(indtop,y(indtop),1);
b2i = pfit(1);                 % slope at upper part of cell
a1i = mean(indcell);           % x position of kink
c1i = mean(y(indcell));   % additive offset

% bounds
lb_b = min([b1i b2i]);
%lb_b = lb_b/2^sign(lb_b);
lb_b = - abs(lb_b)*20;
ub_b = max([b1i b2i]);
%ub_b = ub_b*2^sign(ub_b);
ub_b = abs(ub_b)*20;
start_values = [a1i b1i b2i c1i];

% setting fit options and parameter constraints
fitopt = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[min(indbot) lb_b lb_b min(y(indcell))],...
    'Upper',[max(indtop) ub_b ub_b max(y(indcell))],...
    'Startpoint',start_values);

% fit 
cfun = fit(indcell,y(indcell),fit_model,fitopt);


indkink = round(cfun.a1);
slopetop = cfun.b2;
slopebot = cfun.b1;

end

