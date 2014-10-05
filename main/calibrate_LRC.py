'''
Created on Sep 9, 2014

@author: Yigong
'''
from tables import openFile

base_directory = '/Users/Yigong/Google Drive/Research/CCD/matlab/calibration/'

energy_top = np.array([])
energy_bottom = np.array([])

for i in range(1,3):
    mat_name_tmp = base_directory + 'tracks_%s.mat' % i
    file_tmp = openFile(mat_name_tmp)
    energy_top_tmp = np.squeeze(file_tmp.root.energy_top_tmp[:])
    energy_bottom_tmp = np.squeeze(file_tmp.root.energy_bottom_tmp[:])
    
    energy_top = np.concatenate((energy_top, energy_top_tmp))
    energy_bottom = np.concatenate((energy_bottom, energy_bottom_tmp))
    
    