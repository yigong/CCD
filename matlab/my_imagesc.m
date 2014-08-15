function my_imagesc(data, color_map)

imagesc(data, [-500,500]);
colormap(color_map);
colorbar;
axis image xy
