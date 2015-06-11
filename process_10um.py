import numpy as np
from PyBeamDiag.G4Run import G4Run
from PyBeamDiag.G4Event import G4Event
from scipy.io import loadmat
import pickle

# cd to data dir
fileDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out/10um'
aRunEventList = []
for fileIdx in range(15):
  fileName = fileDir + '/G4W_%s_10um.out' %fileIdx
  aThread = G4Run(fileName)
  aThread.parse()
  aRunEventList.extend(aThread.eventList) 
pickle.dump(aRunEventList, open(fileDir+'/event_list.p', "wb"))

# load psf
psfDir = '/global/scratch/ygzhang/PyBeamDiag/data/psf.mat'
psfTable = []
data_tmp = loadmat(psfDir) 
for i, p in enumerate(data['psf']):
  psfTable.append(p[30:91, 30:91]) # using 60 x 60 psf
