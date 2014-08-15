from CCD.calibrate import self_subtraction
from CCD.io import save_mat

folder = '/Users/Yigong/Research/photon_diagnostics/data/08_12_14_beta_gun/'
    #folder = '/Users/Yigong/Documents/beamCCD/data/07_22_14_200K_Cs137/'
    #folder = "/Volumes/TEAM 7B/07_22_14_200K_Cs137/"
    
file_names = []
for image_index in range(8, 14):
    file_name = folder + 'image_%s.fit' % image_index
    print image_index
    median_subtracted_tmp = self_subtraction(file_name)
    variable_name = 'median_subtracted_%s' % image_index
    file_name_mat = variable_name
    save_mat(variable_name, median_subtracted_tmp, file_name_mat)