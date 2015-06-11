import numpy as np
from PyBeamDiag.G4Run import G4Run
from PyBeamDiag.G4Event import G4Event
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
