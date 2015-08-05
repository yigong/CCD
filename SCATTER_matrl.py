import sys
import os
import multiprocessing
from functools import partial
from PyBeamDiag.g4Parser import parse_E_phi

def main():

    # arguments
    args = sys.argv
    outDir = args[1]

    G4Dir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/mat'
    chdir(G4Dir)
    fileList = glob('./G4W_*_Li.out')
    
    Ephi_partial = partial(parse_E_phi, outFolder=outDir)
    
    pool = multiprocessing.Pool(15)
    pool.map(Ephi_partial, fileList)

if __name__ == "__main__":
    print "start to parse E and phi"
    main()


