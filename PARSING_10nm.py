import os
import glob
import multiprocessing
from PyBeamDiag.g4Parser import parse
from functools import partial

@profile
def main():
    # cd to G4out dir
    dataDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/10nm/split/split'
    os.chdir(dataDir)
    fileList = glob.glob('*.G4track')
    
    ## serial
    #parse(fileList[0])
    
    # parallel
    pool = multiprocessing.Pool(15)
    pool.map(partial(parse, psfTable), fileList)

if __name__ == "__main__":
    main()
    print 'a script is being run'
