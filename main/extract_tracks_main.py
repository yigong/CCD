from CCD.calibrate.calibrate import calculate_median_image
from CCD.extract_tracks import extract_tracks

folder = '/Users/Yigong/Documents/beamCCD/data/08_05_14_background_10min/'

file_names = []
for image_index in range(1,21):
    file_name_temp = folder + 'image_%s.fit' % image_index
    file_names.append(file_name_temp)
    
median_image, median_subtracted_image_stack = calculate_median_image(folder, file_names,
                                                   save_median_flag=True,
                                                   save_median_subtracted_flag=True)

track_list, energy_list, length_list, position_list = extract_tracks(
                                                        median_subtracted_image_stack, 
                                                        5.736e-4, 800)

#track_list, energy_list = init_segment_calibrate(median_subtracted_stack)



