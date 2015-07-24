from scipy.io import savemat
from pickle import load
from os import chdir
from glob import glob
from re import split

dataDir = '/Users/Yigong/Google Drive/Research/CCD/data/step_size'
chdir(dataDir)
fileList = glob('PSF_50um_*P.p')
for fName in fileList:
    preName = split('\.', fName)[0]
    data_list = load(open(fName, 'rb'))
    keys = ['xWindow', 'alphaT', 'alphaM', 'betaM', 'xT', 'yT',
           'xM', 'yM', 'ET', 'EM']
    data_dict = dict([pair for pair in zip(keys, data_list)]) 
    savemat('%s.mat' %(preName), data_dict)


