import os
import glob
import multiprocessing
from PyBeamDiag.G4Parser import parse

# cd to G4out dir
dataDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/3um'
os.chdir(dataDir)
fileList = glob.glob('./G4W_*_3um.out')

# parallel
pool = multiprocessing.Pool(15)
pool.map(parse, fileList)

if __name__ == "__main__":
    print 'a script is being run'
