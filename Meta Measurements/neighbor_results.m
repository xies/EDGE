function out = neighbor_results(data, name)

times = 1:size(data, 1);
bins = 3:9;

n = hist(data.', bins);
n = bsxfun(@rdivide, n, sum(n));  % normalize

%%
figure(1);
for i = times
    bar(bins, n(:,i)');
    xlabel(['# of neighbors (' name ')']);
    ylabel('frequency')
    title(['time = ' num2str(i)]);
    axis([min(bins) max(bins) 0 max(n(:))]);
    pause(0.05);
end
%%
figure(2);
out = my_std(data);
plot(times, out);
xlabel('time')
ylabel(['standard deviation of # of neighbors (' name ')']);
pause(5)

