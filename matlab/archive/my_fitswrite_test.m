%folder_in = '/Users/Yigong/Research/photon_diagnostics/data/08_08_14_background/raw/';
folder_in = '/Users/Yigong/Research/photon_diagnostics/data/08_18_14_Am241/measurement/';
folder_out = '/Users/Yigong/Research/photon_diagnostics/data/08_18_14_Am241/cropped/';
for i = 1:100
    filename_tmp = sprintf('image_%s.fit', int2str(i))
    filename_in = [folder_in, filename_tmp]
    filename_out = [folder_out, filename_tmp]
    my_fitswrite(filename_in, filename_out)
end


