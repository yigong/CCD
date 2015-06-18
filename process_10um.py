from os import chdir
from glob import glob

# initialize variables
endOfFile = False

# cd to data directory
fileDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/10um'
chdir(fileDir)

# find file names
outFilenames = glob("./G4W_*_10um.out")

for outFile in outFilenames:
  aRun = G4Run(outFile)
  while not endOfFile:
    aRun.parse1Event()
    endOfFile = aRun.endFlag
     
