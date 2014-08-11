

# for 70A
# data_dir = '/Users/Yigong/Documents/beamCCD/data/'
# folder = data_dir + '07_21_14_200K_Cali/'

# for MacBook
folder = '/Users/Yigong/Documents/beamCCD/data/07_25_14_Y88/'
# @folder, the folder where the data are stored

fnames_in = []
for img_idx in range(10, 40):
    fileName = folder + 'Y88_image%s.fit' % img_idx
    fnames_in.append(fileName) 
# @var fnames_in, abosolute address of the data 
    
load_flag = True
segment_flag = True
calibrate_flag = False
# @var load_flag, the flag to load the data or not
# @var segment_flag, the flag to segment the tracks or not.
# @var calibrate_flag, the flag to calibrate the tracks or not.

if load_flag == True:
    [medianImage, imageStack] = median_img(folder, fnames_in, saveMedian_flag=True, saveMedianSubtracted_flag=True)
    medianSubtractedImageStack = imageStack - medianImage
if segment_flag == True:    
    track_list, energy_list = segment(medianSubtractedImageStack, n_images=len(fnames_in))
if calibrate_flag == True:
    track_calibrated_list, energy_calibrated_list, gain = calibrate(track_list, energy_list, 60, plot_flag=True)
