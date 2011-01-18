function position = get_nuclear_position(intensity, t, c)

% figure(1)
% for i = 1:size(intensity, 1)
    x = intensity(t, :, c); 
    x=x(:);
%     plot(x); 
%     hold on
    maxi = find( diff([NaN; x])>0 &  diff([x; NaN])<0   );
    % of these "peaks", find the one at the maximum intensity
    
    if length(maxi) > 1
        tosort = [maxi(:) x(maxi)];
        tosort = sortrows(tosort);
        maxi = tosort(end,1);
    end
    
    if ~isempty(maxi)
        % just a sanity check: the peak should be near the top!
        if x(maxi) < my_mean(x')
            maxi = [];
        end
        if sum(~isnan(x)) < 5  % want at least 5 good points to get an estimate
            maxi = [];
        end
    end
    
%     plot(maxi, x(maxi), '.r');
%     hold off
%     pause(1.2);
% end

if isempty(maxi)
    position = NaN;
else
    position = maxi;
end