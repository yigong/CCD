function dummy = update_image(source_name, destination_name, color_map)
folder = '/Users/Yigong/Google Drive/Research/CCD/data/';
data = fitsread([folder, source_name]);

imagesc(data, [0,2^16]);
colormap(color_map);
colorbar;
axis image xy

imwrite(data/2^16*255, color_map, sprintf('./%s',destination_name) , 'png')
