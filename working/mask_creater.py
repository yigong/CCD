import numpy as np

def create_mask(image_dimension):
    '''
    generate the mask to extract effective area for a certain full image dimension
    '''
    n_row = image_dimension[0]
    n_col = image_dimension[1]
    mask = np.zeros(image_dimension)
    mask[2:1755, 12:3518] = 1
    mask[(n_row-2-1752-1):(n_row-2), 12:3518] = 1
    mask.astype(bool)
    return mask
    