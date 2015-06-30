import os
import glob
import multiprocessing
from PyBeamDiag.g4Parser import parse

# cd to G4out dir
dataDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/1um'
os.chdir(dataDir)
fileList = glob.glob('./G4W_*_1um.out')

# parallel
pool = multiprocessing.Pool(15)
pool.map(parse, fileList)

if __name__ == "__main__":
    print 'a script is being run'
