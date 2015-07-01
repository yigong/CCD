import numpy as np
from numpy import tan, deg2rad
from tables import openFile
from glob import glob
from pickle import dump
from os import chdir

dataDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/1um/mat'
chdir(dataDir)
Files = glob('./*.mat')
numTracks = len(Files)

recon = np.zeros(4)
x_CCD = np.zeros(numTracks)
y_CCD = np.zeros_like(x_CCD)
alpha = np.zeros_like(x_CCD)
beta  = np.zeros_like(x_CCD)
x_window = np.zeros_like(x_CCD)
print 'after initialization'

for i, matFile in enumerate(Files):
    file = openFile(matFile)
    recon = np.squeeze(file.root.recon)
    x_CCD[i] = 10.5/2 + recon[2] * 10.5 - 37000./2
    y_CCD[i] = 10.5/2 + recon[3] * 10.5 + 2000.
    alpha[i] = recon[0]
    beta[i]  = recon[1]

print 'finish loading data'
x_window = x_CCD - y_CCD * tan(deg2rad(alpha)) 
x_window[(90<alpha) & (alpha<270)] = float('inf')
dump(x_window, open('../x_window.p', 'wb'))
dump([x_window, x_CCD, y_CCD, alpha, beta], open('../ridgeFollow.p', 'wb'))
print 'data is dumped'
