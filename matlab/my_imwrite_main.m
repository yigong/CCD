 
folder_in = '/Users/Yigong/Research/photon_diagnostics/data/10_27_14_doubleThickness/measurement/';
for i = 1
filename_tmp = sprintf('image_%s_corrected.fit', int2str(i));
filename_in = [folder_in, filename_tmp]
image_tmp = fitsread(filename_in);
sprintf('image_%s', int2str(i))
my_imwrite(image_tmp, cmaphotlog, sprintf('image_%s', int2str(i)))
end