from CCD.calibrate import self_subtraction


image_index = 1
folder = '/Users/Yigong/Research/photon_diagnostics/data/08_18_14_Am241/measurement/'
file_name = folder + 'image_%s.fit' % image_index

self_subtraction(file_name, save_flag=True)

