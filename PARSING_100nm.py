import os
import glob
import multiprocessing
import numpy as np
from scipy.io import loadmat
from PyBeamDiag.g4Parser import parse
from functools import partial

#@profile
def main():
    # get current working dir
    CWD = os.getcwd()
    # cd to G4out dir
    dataDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/100nm/split'
    os.chdir(dataDir)
    fileList = glob.glob('*.G4track')
    
    # load PSF
    psfDir = '/global/scratch/ygzhang/PyBeamDiag/data/psf.mat'
    psfTable = np.array(loadmat(psfDir)['psf'])

    ## serial
    map(partial(parse, psfTable), fileList[:10])
    
    # parallel
    pool = multiprocessing.Pool(15)
    pool.map(partial(parse, psfTable), fileList)

    # go back to initial directory
    os.chdir(CWD)

if __name__ == "__main__":
    main()
    print 'a script is being run'
