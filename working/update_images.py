'''
Created on Aug 8, 2014

@author: Yigong

Task:   
update the images being put in google drive. The upper and lower quadrants
are matched.
'''

from CCD.io import load_fits, save_fits
from CCD.plot import image2d

images = ['Y88.fit', 'Am241.fit', 'Cs137.fit', 'background.fit']

for image_name in images:
    folder = '/Users/Yigong/Desktop/'
    image_tmp = load_fits(folder+image_name, True)
    image_tmp_low = image_tmp[:1753, :]
    image_tmp_up = image_tmp[1761:, ...]
    image = np.concatenate((image_tmp_low, image_tmp_up), axis=0)
    save_fits(image, folder+image_name, True)