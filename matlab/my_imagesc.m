function my_imagesc(data, color_map)

imagesc(data, [0,2^16]);
colormap(color_map);
colorbar;
axis image xy
