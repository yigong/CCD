function my_imwrite(image, color_map, file_name)
imwrite(image/2^16*255, color_map, sprintf('./%s',[file_name,'.png']) , 'png')
