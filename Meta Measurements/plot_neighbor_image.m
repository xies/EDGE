function plot_neighbor_image(data, fignum)

figure(fignum);
img = data';
img(isnan(img)) = -3;
imagesc(img);
colorbar;
xlabel('time');
ylabel('cell index');