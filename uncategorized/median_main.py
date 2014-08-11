data_dir = '/Users/Yigong/Documents/beamCCD/data/'

folder = data_dir + '04_28_14_Am_calibration/'
img_idx = range(1,101)

[med_img, img_array] = median_img(folder=folder, img_idx=img_idx, save_flag=False)