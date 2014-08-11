from io_func import *

def median_img(folder, fileNames_in, save_median_flag, save_median_subtracted_flag):
    """get median image across a set of images.
    Args:
         folder: a string of the folder where the data are stored
         fileNames_in: a list of the absolute address of the images
         save_med_flag: a flag of whether to save the median image
         save_med_sub_flag: a flag of whether to save the median subtracted images
    
    Returns:
         medianImage: an array of median image
         imageStack: an 3D array of median subtracted images.
        
    Raises:
    """
    import astropy.io.fits as pyfits
    import numpy as np
    from datetime import datetime
    
    row_num = 3514
    col_num = 3506
    # define the size of the array
        
    imageStack = np.zeros((len(fileNames_in), row_num, col_num), dtype=np.uint16)
    # create a 3d image stack
    
    # what data type to use?
    # uint16: unsigned integer (0 to 65536)
    # int16 : integer (-32768 to 32767)
    
    for idx, fileName in enumerate(fileNames_in):
        
        full_img = load_fits(fileName)
        # load the image data into @var full_img
        
        img_data = np.array(full_img[2:3516, 12:3518], dtype=np.uint16)
        # fill the image stack
        
        imageStack[idx, ...] = img_data
        img_data = np.zeros_like(img_data)
        # clear the temp image
        
        print 'finish loading image # %s out of %s into stack' %(idx, len(fileNames_in))
        # print the progress
        #1/0

    medianImage = np.median(imageStack, 0)
    # take the median 
    
    if save_median_flag == True:
        # save the median image
        median_name = 
        saveName = folder+'median%s.fit'
        save_fits(saveName, medianImage)
         
    if save_median_subtracted_flag == True:
        # save median subtracted images
        for idx, fileName in enumerate(fileNames_in):
            suffix = fileName.split('image')[-1]
            saveName = folder+'medianSubtracted_%s' % suffix
            imageTemp = imageStack[idx, ...] - medianImage
            save_fits(saveName, imageTemp)
    
    return [medianImage, imageStack]
    
def medianEZ_img(folder, img_idx, plot_flag, save_flag):
    """ the median value of a pixel is the median of that row.
    
    Args:
        folder
        img_idx
        plot_flag
        save_flag
    Returns:
        asdfasdf
    Raises:
    """
    
    import astropy.io.fits as pyfits
    import numpy as np
    
    fname = folder + 'image_' + str(img_idx) + '.fit'
    data_cube = pyfits.getdata(fname, 0)
    img_U1 = data_cube[1766:, 1773:3524]
    med_img = np.zeros_like(img_U1)
    nrow, ncol = np.shape(img_U1)
    for row_idx in range(nrow):
        med_img[row_idx, :] = np.median(img_U1[row_idx, :])
    return [med_img, img_U1]
    if save_flag == True:
            hdu = pyfits.PrimaryHDU(med_img)
            hdulist = pyfits.HDUList([hdu])
            hdulist.writeto(folder+'median_ez.fit')    
    return (med_img, imageStack)
    
        
    
        
        
        