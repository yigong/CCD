import numpy as np
import matplotlib.pyplot as plt
from pylab import *

from scipy.optimize import curve_fit
from gaussianFit_func import gaussian_noBG


# import fits related modules and curve fitting function
import astropy.io.fits as pyfits

def histogram_img(folder, img_idx, plot_flag):
    baseline = list()
    noise = list()
    
    for i in img_idx:
        if i < 10:
            fname = folder+'image_' + str(i)+'.fit'
        elif  10 <= i < 100:
            fname = folder+'image_' + str(i)+'.fit'
        else:
            fname = folder+'image_' + str(i)+'.fit'
        # data extraction
        data_cube = pyfits.getdata(fname, 0)
        img = data_cube[1766:, 1773:3524]
        BN = 200
        RG = (0, 2000)
        if plot_flag:
            fig1, [ax1, ax2] = plt.subplots(1,2)
            ax1.imshow(img, cmap=plt.cm.hot, interpolation='nearest')
            ax1.set_title('image_'+str(i))
            img_1d = img.flatten()
            #1/0
            img_counts, img_bins_temp = np.histogram(img_1d, bins=BN, range=RG)
            img_bins = img_bins_temp[1:] - (img_bins_temp[1]-img_bins_temp[0])/2
        # fitting
        guess = [size(img_1d), median(img_1d), (median(img_1d)-min(img_1d))/3]
        parameter, parameter_cov = curve_fit(gaussian_noBG, img_bins, img_counts, p0=guess)
        if plot_flag:
            ax2.hist(img_1d, bins=BN, range=RG)
            ax2.plot(img_bins, gaussian_noBG(img_bins, *parameter), 'rx-')
#             ax2.set_yscale('log')
            ax2.set_title('baseline='+str(parameter[1])+' || '+'std='+str(parameter[2]))
            plt.show()
        baseline.extend([parameter[1]])
        noise.extend([parameter[2]])
    return [baseline, noise]

def histogram_single_img(fname):
    data_full_img = pyfits.getdata(fname, 0)
    row_num = 3508
    col_num = 3512
    
    data_L1 = data_full_img[row_num/2:, :col_num/2-1]
    data_L2 = data_full_img[:row_num/2-1, :col_num/2-1]
    data_U1 = data_full_img[row_num/2:, col_num/2:]
    data_U2 = data_full_img[:row_num/2-1, col_num/2:]
    data_list = [data_L1, data_L2, data_U1, data_U2]
    
    fig, ax = plt.subplots(2,2)
    for i, a in enumerate(ax.flatten('F')):
        data_flattened = data_list[i].flatten()
        MIN = np.min(data_flattened)
        MAX = np.max(data_flattened)
        bin_plot = np.arange(MIN, MAX+1)
        bin_hist = np.arange(MIN-1, MAX+1) + 0.5
        counts, temp = np.histogram(data_flattened, bin_hist)
        a.plot(bin_plot, counts)
        a.set_title('median = ' + str(np.median(data_flattened)))
    plt.show()
    return 0
        

def histogram_row(folder, img_idx, row_idx, plot_flag):
    baseline = list()
    noise = list()
    
    #   if img_idx < 10:
    fname = folder+'image_' + str(img_idx[0])+'.fit'
    #   elif  10 <= i < 100:
    #       fname = folder+'image_' + str(i)+'.fit'
    #   else:
    #       fname = folder+'image_' + str(i)+'.fit'
    
    # data extraction
    data_cube = pyfits.getdata(fname, 0)
    img = data_cube[1766:, 1773:3524]
    #   1/0
    BN = 100
    RG = (0, 1000)
    for i in row_idx: 
        img_1d = img[i]
        if plot_flag:
            fig1, axs = plt.subplots(1,2)
            axs[0].plot(img_1d)
            axs[0].set_title('row_'+str(i))
        img_counts, img_bins_temp = np.histogram(img_1d, bins=BN, range=RG)
        img_bins = img_bins_temp[1:] - (img_bins_temp[1]-img_bins_temp[0])/2
        # fitting
        guess = [size(img_1d), median(img_1d), (median(img_1d)-min(img_1d))/3]
        parameter, paremeter_cov = curve_fit(gaussian_noBG, img_bins, img_counts, p0=guess)
        if plot_flag:
            axs[1].hist(img_1d, bins=BN, range=RG)
            axs[1].plot(img_bins, gaussian_noBG(img_bins, *parameter), 'rx-')
            #ax2.set_yscale('log')
            axs[1].set_title('baseline='+str(parameter[1])+' || '+'std='+str(parameter[2]))
            plt.show()
        baseline.extend([parameter[1]])
        noise.extend([parameter[2]])
        print i
    return [baseline, noise]


  

 
      


