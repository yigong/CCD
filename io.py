import pickle
import astropy.io.fits as pyfits
import subprocess
import numpy as np

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
    '''
    Args:
        variable_name: string. The variable name in MATLAB workspace.
        value: variable. The variable in python workspace.
        file_name: string. Don't need .mat
    '''
    from scipy.io import savemat
    savemat(file_name, {variable_name:value})
    
def array_to_file(file_name, data):
    file_obj = open(file_name, 'w')
#     data_dict = { line : value for (line, value) in enumerate(data)}
#     string = str()
#     for i in range(len(data)):
#         tmp = '(%s)%%s\n' % i
#         string += tmp
    data_tuple = tuple(data)
    each_line = '%s\n'
    all_lines = len(data) * each_line
    
    file_obj.write(all_lines % data_tuple)
    return (all_lines % data_tuple)

def read_file(file_name):
    import re
    file_obj= open(file_name, 'r')
    col1 = list()
    col2 = list()
    for line in file_obj:
        [col1_tmp, col2_tmp, dum] = re.split(' ', line) 
        
        col1.append(float(col1_tmp))
        col2.append(float(col2_tmp))
    
    return np.array(col1), np.array(col2)




    
