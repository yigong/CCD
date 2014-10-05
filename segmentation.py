import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from skimage.morphology import binary_dilation
from skimage.measure import label
from skimage.segmentation import clear_border

from CCD.gaussian_fit import gaussian_no_background

def noise_sigma(image, bin_range=(-500, 500), bin_step=10, plot_flag=False):
    
    bins_tmp = np.arange(bin_range[0], bin_range[1]+bin_step, bin_step)
    
    pixel_val = image.flatten();
    counts, bins = np.histogram(pixel_val, bins=bins_tmp)
    bins = bins[1:] - (bins[1] - bins[0])/2
    
    MEDIAN = np.median(pixel_val)
    NUM_PIXEL = np.size(pixel_val)
    MIN = np.min(pixel_val)
    
    fit_guess = [NUM_PIXEL, MEDIAN, (MEDIAN-MIN)/3]
    sigma = 1./counts
    
    fitting_parameter, fitting_parameter_cov = curve_fit(gaussian_no_background,
                                                    bins, counts, p0=fit_guess, sigma=sigma)
    
    centroid = fitting_parameter[1]
    std = fitting_parameter[2]
    
    if plot_flag == True:
        fig, ax1 = plt.subplots()
        ax1.hist(pixel_val, bins=bins_tmp)
        ax1.plot(bins, gaussian_no_background(bins, *fitting_parameter), 'rx-')
        ax1.set_xlabel('track energy (ADC)')
        ax1.set_ylabel('counts')
        ax1.set_title('centroid=%.1f \n std=%.2f' % (centroid, std))
        # histogram on ADC unit
        
#         energy_kev = pixel_val * gain
#         bins_kev = bins_temp * gain
#         ax2.hist(energy_kev, bins=bins_kev)
#         ax2.set_xlabel('track energy (keV)')
#         ax2.set_ylabel('counts')
#         ax2.set_title('Am_241 spectrum \n std=%.2f' % std*gain)
#         # histogram on keV
        plt.show()
    
    result = {}
    result['centroid'] = centroid
    result['STD'] = std


    print 'The histogram is centered at %(centroid).4e with STD = %(STD).4e' % result
    
    return result  
    

def segment_tracks(image, thresholds):
    track = []
    
    NEIGHBOR_ARRAY = np.ones((3, 3))
    
    dim0, dim1 = np.shape(image)
      
    image_bottom = image[:dim0/2, :]
    image_top = image[dim0/2:, :]
      
    image_binary_bottom = image_bottom > thresholds[0]
    image_binary_top = image_top > thresholds[1]
    
    image_binary = np.concatenate((image_binary_bottom, image_binary_top))
    image_binary = clear_border(image_binary)
    image_dilated = binary_dilation(image_binary, NEIGHBOR_ARRAY)
    image_segmented, max_label = label(image_dilated, neighbors=8, background=False,
                                       return_num=True)
    print '%s tracks are segmented.' % (max_label)
    for track_idx in range(1,max_label):
        
        if track_idx%100 == 0:
            print 'Have segmented %s tracks' % track_idx
        rows_in_image, cols_in_image = np.where(image_segmented==track_idx)
        
        row_min = np.min(rows_in_image)
        row_max = np.max(rows_in_image)
        col_min = np.min(cols_in_image)
        col_max = np.max(cols_in_image)
        
        rows_in_track = rows_in_image - row_min
        cols_in_track = cols_in_image - col_min
        
        track_dict = {}
        track_dict['track'] = np.zeros([row_max-row_min+1, col_max-col_min+1])
        track_dict['track'][rows_in_track, cols_in_track] = image[rows_in_image, cols_in_image]
        track_dict['track_1d'] = image[image_segmented == track_idx]
        track_dict['energy'] = np.sum(track_dict['track_1d'])
        track_dict['position'] = (rows_in_image, cols_in_image)
        
        track.append(track_dict)

    return track    
        
         
    
    
    