import os
import glob
import multiprocessing
from PyBeamDiag.g4Parser import parse

@profile
def main():
    # cd to G4out dir
    dataDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/30nm'
    os.chdir(dataDir)
    fileList = glob.glob('./G4W_*_0p03um.out')
    
    ## serial
    #parse(fileList[0])
    
    # parallel
    pool = multiprocessing.Pool(15)
    pool.map(parse, fileList)

if __name__ == "__main__":
    main()
    print 'a script is being run'
