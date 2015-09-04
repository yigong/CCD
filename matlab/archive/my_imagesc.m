function my_imagesc(data, color_map)

imagesc(data, [0,2^16*5.5e-4]);
colormap(color_map);
colorbar;
axis image xy
% axis image: square pixel.
% axis xy:    (0, 0) is at the lower left corner. 
