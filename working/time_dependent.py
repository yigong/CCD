from CCD.calibrate import calculate_median_image
from CCD.calibrate import init_segment_calibrate

load_flag = 1
if load_flag:
    folder = '/Users/Yigong/Research/photon_diagnostics/data/08_06_14_140K/'
    #folder = '/Users/Yigong/Documents/beamCCD/data/07_22_14_200K_Cs137/'
    #folder = "/Volumes/TEAM 7B/07_22_14_200K_Cs137/"
    
    file_names = []
    for image_index in range(550, 594):
        file_name_temp = folder + 'image_%s.fit' % image_index
        file_names.append(file_name_temp)
    
    median_image, median_subtracted_image_stack, image_stack = calculate_median_image(folder, file_names,
                                                        save_median_flag=False,
                                                        save_median_subtracted_flag=False)

    
# figure()
# 
# for i in range(0,3000, 10):
#     plot(median_subtracted_image_stack[:, 1, i] )