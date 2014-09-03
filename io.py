import pickle
import astropy.io.fits as pyfits
import subprocess

def load_fits(file_name, ds9_flag=False, header_out_flag=False):
    ''' load fit file
    '''
    
    full_image = pyfits.getdata(file_name, 0)
    
    if ds9_flag == True:
        subprocess.call(['ds9',file_name])
    if header_out_flag == True:
        return full_image, header_out_flag
    else:
        return full_image

def save_fits(image, fileName, ds9_flag=False, header_in=''):
    '''write 2d array into fit file
    Args:
        fileName: absolute path where to write the image
        image: 2d array
        ds9_flag = flag of showing the image with ds9        
    Returns:
        nothing
    
    ''' 
    if header_in == '':      
        hdu = pyfits.PrimaryHDU(image)
    # else add header into the new fits file
    hdu.writeto(fileName, clobber=True)
    print 'An image has been written to %s' %fileName
       
    if ds9_flag == True:
        subprocess.call(['ds9',fileName]) 
        
    return


def save_pickle(data, file_name):
    pickle.dump(data, open(file_name, 'wb'))
    return

def load_pickle(file_name):
    return pickle.load(open(file_name, 'rb'))

def save_mat(variable_name, value, file_name):
    from scipy.io import savemat
    savemat(file_name, {variable_name:value})