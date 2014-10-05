from CCD.calibrate import self_subtraction

folder = '/Users/Yigong/Research/photon_diagnostics/data/08_18_14_Am241/measurement/'

for image_index in np.arange(1,31):
    file_name = folder + 'image_%s.fit' % image_index

    self_subtraction(file_name, save_flag=True)

