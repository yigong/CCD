from CCD.calibrate import self_subtraction, init_segment_calibrate
from CCD.io import save_mat

folder = '/Users/Yigong/Research/photon_diagnostics/data/08_06_14_140K/'
file_names = []
for image_index in range(550, 594):
    file_names.append(folder + 'image_%s.fit' % image_index)

    #folder = '/Users/Yigong/Documents/beamCCD/data/07_22_14_200K_Cs137/'
    #folder = "/Volumes/TEAM 7B/07_22_14_200K_Cs137/"
median_subtracted_image_stack = np.zeros((len(file_names), 3506, 3506))

for i, file_name in enumerate(file_names):
    print i
    median_subtracted_tmp = self_subtraction(file_name)
    median_subtracted_image_stack[i, ...] = median_subtracted_tmp
#     variable_name = 'median_subtracted_%s' % image_index
#     file_name_mat = variable_name
#     save_mat(variable_name, median_subtracted_tmp, file_name_mat)
init_segment_calibrate(median_subtracted_image_stack)
