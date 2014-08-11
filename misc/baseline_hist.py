import numpy as np
import matplotlib.pyplot as plt
from pylab import *

from scipy.optimize import curve_fit
from gaussianFit_func import gaussian_noBG


# import fits related modules and curve fitting function
import astropy.io.fits as pyfits

data_dir = 'home/photon_diag/Documents/images/'
folder = data_dir + '07_14_14_200K_Calibration/'
img_idx = range(3,55,3)

for i in range(0, len(img_idx), 4):
    fig, ax = plt.subplots(2,2)
    for j, a in enumerate(ax.flatten('C')):
        fname = folder+'image' + str(img_idx[i+j])+'.fit'
        print fname
        data_full_img = pyfits.getdata(fname, 0)
        
        row_num = 3508
        col_num = 3512
        
        data_U2 = data_full_img[:row_num/2-1, col_num/2:]
        data_flattened = data_U2.flatten()
        
        MIN = np.min(data_flattened)
        MAX = np.max(data_flattened)
        
        bin_plot = np.arange(MIN, MAX+1)
        bin_hist = np.arange(MIN-1, MAX+1) + 0.5
        
        counts, temp = np.histogram(data_flattened, bin_hist)
        
        a.plot(bin_plot, counts)
        a.set_title('median = ' + str(np.median(data_flattened)))

plt.show()