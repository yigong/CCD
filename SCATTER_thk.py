import sys
import multiprocessing
from functools import partial
from PyBeamDiag.g4Parser import parse_E_phi
from os import chdir
from glob import glob

def main():

    # arguments
    args = sys.argv
    outDir = args[1]
    t = args[2]
    print outDir

    G4Dir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/thk'
    chdir(G4Dir)
    fileList = glob('G4W_*_thk%sum.out'%(t))
    
    Ephi_partial = partial(parse_E_phi, outFolder=outDir, thk=t)
    
    pool = multiprocessing.Pool(15)
    pool.map(Ephi_partial, fileList)

if __name__ == "__main__":
    print "start to parse E and phi"
    main()


