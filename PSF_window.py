import numpy as np
from numpy import tan, deg2rad
from scipy.io import loadmat
from tables import openFile
from glob import glob
from pickle import dump
from os import chdir, getcwd
import sys

def main():

    # parse input arguments
    args = sys.argv
    stepLen = args[1]       # step length,      eg: 3um
    fitsFdr = args[2]       # fits file folder, eg: fits_backP
    outFile = args[3]       # pickle file name, eg: PSF_window_3um 

    currentDir = getcwd()
    dataDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/%s/%s' %(stepLen, fitsFdr)
    chdir(dataDir)
    Files = glob('*.mat')
    numTracks = len(Files)
    
    alphaT = np.zeros(numTracks)
    alphaM = np.zeros_like(alphaT)
    betaM  = np.zeros_like(alphaT)
    xT     = np.zeros_like(alphaT)
    yT     = np.zeros_like(alphaT)
    xM     = np.zeros_like(alphaT)
    yM     = np.zeros_like(alphaT)
    ET     = np.zeros_like(alphaT)
    EM     = np.zeros_like(alphaT)
    
    print 'after initialization'
    for i, matFile in enumerate(Files):
        result = loadmat(matFile, squeeze_me=True)
        alphaT[i] = result['alphaT']
        alphaM[i] = result['alphaM']
        betaM[i] = result['betaM']
        xT[i] = result['xT']
        yT[i] = result['yT']
        xM[i] = result['xM']
        yM[i] = result['yM']
        ET[i] = result['ET']
        EM[i] = result['EM']
    print 'finish loading data'
    xWindow = xM - yM * tan(deg2rad(alphaM)) 
    xWindow[(90<alphaM) & (alphaM<270)] = float('inf')
    dump([xWindow, alphaT, alphaM, betaM, xT, yT, xM, yM, ET, EM], open('../%s.p' %(outFile), 'wb'))
    print 'data is dumped'
    chdir(currentDir)

if __name__ == '__main__':
    print 'IN PSF_window.py'
    main()
