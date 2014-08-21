from CCD.calibrate import calculate_median_image
from CCD.extract_tracks import extract_tracks


folder = '/Users/Yigong/Research/photon_diagnostics/data/08_18_14_Am241/measurement/'
#folder = '/Users/Yigong/Documents/beamCCD/data/07_22_14_200K_Cs137/'
#folder = "/Volumes/TEAM 7B/07_22_14_200K_Cs137/"

file_names = []
for image_index in range(1, 101):
    file_name_temp = folder + 'image_%s.fit' % image_index
    file_names.append(file_name_temp)
    
median_image, median_subtracted_image_stack, image_stack = calculate_median_image(folder, file_names,
                                                   save_median_flag=True,
                                                   save_median_subtracted_flag=False)

image_stack_low = median_subtracted_image_stack[:,0:1753,:]
image_stack_up = median_subtracted_image_stack[:,1753:3506,:]

track_list, energy_list_low, length_list_low, position_list = extract_tracks(
                                                        image_stack_low, 1, 500)

track_list, energy_list_up, length_list_up, position_list = extract_tracks(
                                                        image_stack_up, 1, 500)

                                                                     
