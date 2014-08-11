import numpy as np
import matplotlib.pyplot as plt

from skimage.morphology import label, binary_dilation
import astropy.io.fits as pyfits

def segment(image_stack,  threshold=500):
    track_list = []
    energy_list = []
    neighbour_array = np.ones((3,3))

    for i in range(n_images):
        print 'calibrate image #'+ str(i) + ' out of' + str(n_images)
        # find the median subtracted image
        med_sub = image_stack[i, ...]
        # get the binary image for segmentation
        img_binary = med_sub > threshold
        # build the matrix for dialation
        img_dialated = binary_dilation(img_binary, neighbour_array)
        # segmentation
        img_segmented = label(img_dialated, 8, 0)
        # get the total number of tracks
        n_tracks = np.max(img_segmented)
        1/0
        for track_idx in np.arange(n_tracks):
            # get the region of the ith track
            rows, cols = np.where(img_segmented==track_idx)
            rows_min = np.min(rows)
            rows_max = np.max(rows)
            cols_min = np.min(cols)
            cols_max = np.max(cols)
            # crop the track out
            track_temp = med_sub[rows_min:rows_max, cols_min:cols_max]
            energy_list.append(np.sum(track_temp))
            track_list.append(track_temp)
            track_temp = []
            
    return track_list, energy_list

def calibrate(track_list, energy_list, cali_energy, plot_flag=True):
    from scipy.optimize import curve_fit
    from gaussianFit_func import *
    
    # setup bin range and bin width
    RG = (53600, 57700)
    step_adu = 100
    BN_adu = np.arange(RG[0], RG[1], step_adu)
    # histogram
    counts, bins_temp = np.histogram(energy_list, bins=BN_adu)
    bins_adj = bins_temp[1:] - (bins_temp[1]-bins_temp[0])/2
    energy_array = np.array(energy_list)
    # crop the histogram based on the range, RG.
    in_range_bool = np.logical_and(energy_array>RG[0],energy_array<RG[1])
    energy_in_range = energy_array[in_range_bool]
    guess = [np.size(energy_in_range), np.median(energy_in_range), (np.median(energy_in_range)-np.min(energy_in_range))/3]

    sigma = 1./counts
    parameter, parameter_cov = curve_fit(gaussian_noBG, bins_adj, counts, p0=guess, sigma=sigma)
    if plot_flag == True:
        f2, ax2 = plt.subplots()
        ax2.hist(energy_list, bins=BN_adu, range=RG)
        ax2.plot(bins_adj, gaussian_noBG(bins_adj, *parameter), 'rx-')
        #             ax2.set_yscale('log')
        ax2.set_title('baseline='+str(parameter[1])+' || '+'std='+str(parameter[2]))
        ax2.set_yscale('log')
        plt.show()
    gain = cali_energy/parameter[1]
    track_calibrated_list = []
    energy_calibrated_list = []
        
    for i in range(len(energy_list)):
        if i % 1000 == 0:
            print i
        
        track_calibrated_list.append(track_list[i]*gain)
        energy_calibrated_list.append(energy_list[i]*gain)
    
    if plot_flag == True:
        f1, ax1 = plt.subplots()
        BN_keV = np.arange(0, RG[1]*gain, step_adu*gain)
        ax1.hist(energy_calibrated_list, bins=BN_keV)
        ax1.set_title('Am_241 Spectrum || '+'std='+str(parameter[2]*gain))
        plt.show()
   


    return track_calibrated_list, energy_calibrated_list, gain

