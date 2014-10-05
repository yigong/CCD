 
folder_in = '/Users/Yigong/Research/photon_diagnostics/data/09_02_14_Co60/measurement/';
for i = 1:10
filename_tmp = sprintf('image_%s_corrected.fit', int2str(i));
filename_in = [folder_in, filename_tmp];
image_tmp = fitsread(filename_in);
my_imwrite(image_tmp, cmaphotlog, sprintf('image_%s', int2str(i)))
end