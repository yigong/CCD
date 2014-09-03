function my_fitswrite(filename_in, filename_out)
data = fitsread(filename_in);
[row_number, col_number] = size(data)
low_tmp = data(3:1755, 13:3518);
up_tmp = data(row_number-2-1752 : row_number-2, 13:3518);
image = cat(1, low_tmp, up_tmp);
image = int32(image);
fitswrite(image, filename_out)


    



end
