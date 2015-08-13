import os
import glob
import multiprocessing
import numpy as np
from scipy.io import loadmat
from PyBeamDiag.g4Parser import parse
from functools import partial
import sys

#@profile
def main():

    # arguments
    args = sys.argv
    stepLen = args[1]         # step length,            eg: 3um
    fitsFdr = args[2]         # fits file folder,       eg: fits_backP 
    pfacing = args[3]         # pixel-plane facing,     eg: top | bottom 

    # get current working dir
    CWD = os.getcwd()
    # cd to G4out dir
    splitDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/%s/split' %(stepLen)
    os.chdir(splitDir)
    fileList = glob.glob('*.G4track')
    
    # load PSF
    psfDir = '/global/scratch/ygzhang/PyBeamDiag/data/psf.mat'
    psf = np.array(loadmat(psfDir)['psf'])

    ## serial
    ## map(partial(parse, psfTable), fileList[:10])
    
    # partial func
    parse_partial = partial(parse, psfTable=psf, pixelPlane=pfacing, outFolder='%s/%s' %(stepLen, fitsFdr))

    # parallel
    pool = multiprocessing.Pool(15)
    pool.map(parse_partial, fileList)

    # go back to initial directory
    os.chdir(CWD)

if __name__ == "__main__":
    print 'IN PARSING.py'
    main()
