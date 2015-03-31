from CCD.calibrate import self_subtraction

folder = '/Users/Yigong/Research/photon_diagnostics/data/10_27_14_doubleThickness/measurement/'

for image_index in np.arange(1,2):
    file_name = folder + 'image_%s.fit' % image_index

    a = self_subtraction(file_name, save_flag=True)

