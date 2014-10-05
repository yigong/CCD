import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

from CCD.gaussian_fit import gaussian_no_background

from CCD.io import load_fits
from CCD.io import save_fits


def calculate_median_image(save_median_dir, file_names, save_median_flag=True, save_median_subtracted_flag=False):
    """calculate median image across a set of images.
    Args:
         save_median_dir: the folder to save the median image
         file_names: a list of the absolute address of the images
         save_med_flag: a flag of whether to save the median image
         
    Returns:
         median_image: an array of median image
        
    Raises:
    """
    # define the size of the array
    effective_row_number = 3506
    effective_col_number = 3506
    
    image_stack = np.zeros((len(file_names), effective_row_number, effective_col_number), 
                           dtype=np.int32)  

    # create a 3d image stack
    
    
    # what data type to use?
    # uint16: unsigned integer (0 to 65536)
    # int16 : integer (-32768 to 32767)
    
    for i, file_name in enumerate(file_names):
        
        full_img = load_fits(file_name)
        if i == 0:
            row_number, col_number = np.shape(full_img)
        # load the image data into @var full_img
        low_tmp = full_img[2:1755, 12:3518]
        up_tmp = full_img[(row_number-2-1752-1):(row_number-2), 12:3518]
        image_tmp = np.concatenate((low_tmp, up_tmp), axis=0)
    
        image_stack[i, ...] = image_tmp
        # fill the image stack
                
        image_tmp = np.zeros_like(image_tmp)
        # clear the temp image
        
        print 'finish loading image # %s out of %s into stack' %(i+1, len(file_names))
        # print the progress
        #1/0

    median_image = np.median(image_stack, 0).astype(np.int32)
    # take the median 
    
    if save_median_flag == True:
        # save the median image
        median_suffix = file_names[0].split('image')[-1]
        save_name = save_median_dir+'median%s' %median_suffix
        save_fits(median_image, save_name)
    
    

    if save_median_subtracted_flag == True:
        # save the median subtracted image
        median_subtracted_image_stack = image_stack - median_image

        for idx, name in enumerate(file_names):
            suffix = name.split('image')[-1]
            save_name = save_median_dir + 'median_subtracted%s' % suffix
            image_tmp = median_subtracted_image_stack[idx, ...]
            
            save_fits(image_tmp, save_name)
    else:
        median_subtracted_image_stack = image_stack - median_image
    
    return (median_image, median_subtracted_image_stack, image_stack)   

def init_segment_calibrate(image_stack, energy_list=[], threshold=500, plot_flag=True):
    '''segment the tracks for calibration
    
    Args:
    
        image_stack: ~
        threshold: threshold for initial segmentation, in ADC value
        
    Return:
        calibrate_coeff: coefficient to multiply in calibration

    '''
    from skimage.morphology import binary_dilation
    from skimage.measure import label
    
    
    calibrate_energy = 59.5
    
    if energy_list == []:
        track_list = []
        energy_list = []
        neighbour_array = np.ones((3,3))
        
        # create track_list and energy_list to store the tracks and its total energy
        
        for image_idx in range(len(image_stack)):
            image_tmp = image_stack[image_idx, ...]
            image_binary = image_tmp > threshold
            # get the binary image with the givin threshold
            image_dilated = binary_dilation(image_binary, neighbour_array)
            # dilate the image 
            image_segmented, max_label = label(image_dilated, neighbors=8, 
                                               background=-1, return_num=True)
            # segment the tracks and label them with number
            print '%s tracks in image # %s(out of %s)' % (max_label, image_idx+1, len(image_stack)) 
                     
            for track_idx in range(1,max_label+1):
                #print '  NO.%s out of %s tracks' %(track_idx, max_label)
                
                track_temp = image_tmp[image_segmented == track_idx]
                
                energy_temp = np.sum(track_temp)
                #1/0
                energy_list.append(energy_temp)
    1/0    

        
    bin_range = (100000, 107000)
    bin_step = 100
    bins_temp = np.arange(bin_range[0], bin_range[1]+bin_step, bin_step)
    # setup the bins
    
    counts, bins = np.histogram(energy_list, bins=bins_temp)
    # histogram
    
    bins = bins[1:] - (bins[1] - bins[0])/2
    energy = np.array(energy_list)
    
    in_range_bool = np.logical_and(energy>bin_range[0], energy<bin_range[1])
    energy_in_range = energy[in_range_bool]
    # crop the histogram based on the range
    
    fit_guess = [np.size(energy_in_range), np.median(energy_in_range), 
                 (np.median(energy_in_range)-np.min(energy_in_range))/3]
    sigma = 1./counts
    fitting_parameter, fitting_parameter_cov = curve_fit(gaussian_no_background,
                                                        bins, counts, p0=fit_guess, sigma=sigma)
    centroid = fitting_parameter[1]
    std = fitting_parameter[2]
    # fit the histogram with guassian curve
    
    gain = calibrate_energy/centroid
    
    
    if plot_flag == True:
        fig, (ax1, ax2) = plt.subplots(1,2)
        ax1.hist(energy, bins=bins_temp)
        ax1.plot(bins, gaussian_no_background(bins, *fitting_parameter), 'rx-')
        ax1.set_xlabel('track energy (ADC)')
        ax1.set_ylabel('counts')
        ax1.set_title('centroid=%.1f \n std=%.2f' % (centroid, std))
        # histogram on ADC unit
        
        energy_kev = energy * gain
        bins_kev = bins_temp * gain
        ax2.hist(energy_kev, bins=bins_kev)
        ax2.set_xlabel('track energy (keV)')
        ax2.set_ylabel('counts')
        ax2.set_title('Am_241 spectrum \n std=%.2f' % std*gain)
        # histogram on keV
        plt.show()

    return gain
            
def self_subtraction(file_name, save_flag=False, suffix='_corrected'):
    full_img = load_fits(file_name)
    row_number, col_number = np.shape(full_img)
    # load the image data into @var full_img
    low_tmp = full_img[2:1755, 12:3518]
    low_row_median_subtracted, low_row_median = row_median_subtraction(low_tmp)
    low_col_median_subtracted, low_col_median = col_median_subtraction(low_row_median_subtracted) 
    
    up_tmp = full_img[(row_number-2-1752-1):(row_number-2), 12:3518]
    up_row_median_subtracted, up_row_median = row_median_subtraction(up_tmp)
    up_col_median_subtracted, up_col_median = col_median_subtraction(up_row_median_subtracted) 
    
    median_subtracted = np.concatenate((low_col_median_subtracted, up_col_median_subtracted), axis=0)
    
    if save_flag == True:
        # add suffix to the save_file_name
        split_list = file_name.split('.')
        save_file_name = split_list[0] + suffix + '.fit'
        save_fits(median_subtracted, save_file_name, ds9_flag=False)
        
    return median_subtracted
        
def row_median_subtraction(image):
    '''return the row median subtracted image
    '''
    row_median = np.median(image, axis=1)
    row_median_subtracted = np.zeros_like(image)
    
    for i in range(len(image)):
        row_median_subtracted[i, :] = image[i, :] - row_median[i]
        
    return row_median_subtracted, row_median 

def col_median_subtraction(image):
    
    col_median = np.median(image, axis=0)
    col_median_subtracted = np.zeros_like(image)
    
    for i in range(len(image[0])):
        col_median_subtracted[:, i] = image[:, i] - col_median[i]
        
    return col_median_subtracted, col_median

def Am241_calibration(energy, plot_flag=False):
    CALIBRATE_ENERGY = 59.54
    # energy is an array
    bin_range = (104000, 114000)
    bin_step = 100
    bins_temp = np.arange(bin_range[0], bin_range[1]+bin_step, bin_step)
    # setup the bins
    
    counts, bins = np.histogram(energy, bins=bins_temp)
    # histogram
    
    bins = bins[1:] - (bins[1] - bins[0])/2
    
    
    in_range_bool = np.logical_and(energy>bin_range[0], energy<bin_range[1])
    energy_in_range = energy[in_range_bool]
    # crop the histogram based on the range
    
    fit_guess = [np.size(energy_in_range), np.median(energy_in_range), 
                 (np.median(energy_in_range)-np.min(energy_in_range))/3]
    sigma = 1./counts
    fitting_parameter, fitting_parameter_cov = curve_fit(gaussian_no_background,
                                                        bins, counts, p0=fit_guess, sigma=sigma)
    centroid = fitting_parameter[1]
    std = fitting_parameter[2]
    FWHM = std * 2.35
    # fit the histogram with guassian curve
    
    gain = CALIBRATE_ENERGY/centroid
    
    
    if plot_flag == True:
        fig, (ax1, ax2) = plt.subplots(1,2)
        ax1.hist(energy, bins=bins_temp)
        ax1.plot(bins, gaussian_no_background(bins, *fitting_parameter), 'rx-')
        ax1.set_xlabel('track energy (ADC)')
        ax1.set_ylabel('counts')
        ax1.set_title('centroid=%.1f \n FWHM=%.2f' % (centroid, FWHM))
        # histogram on ADC unit
        
        energy_kev = energy * gain
        bins_kev = bins_temp * gain
        ax2.hist(energy_kev, bins=bins_kev)
        ax2.set_title('std=%.2f' % std*gain)

        ax2.set_xlabel('track energy (keV)')
        ax2.set_ylabel('counts')
        # histogram on keV
        plt.show()
    
    result = {}
    result['centroid'] = centroid
    result['FWHM_ADC'] = FWHM
    result['gain'] = gain
    result['FWHM_keV'] = FWHM * gain

    print 'The histogram is centered at %(centroid).4e with FWHM = %(FWHM_ADC).4e' % result
    print 'The CCD has a gain of %(gain).4e' % result
    print 'FWHM at 59.5 keV is %(FWHM_keV).4e keV' % result 

    return result










    
    