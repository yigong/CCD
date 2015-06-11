import numpy as np
from PyBeamDiag.G4Run import G4Run
from PyBeamDiag.G4Event import G4Event
from PyBeamDiag.io import save_fits
from scipy.io import loadmat
import pickle

### cd to data dir
fileDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out/10um_test'
aRunEventList = []
fileIdx = 3
fileName = fileDir + '/G4W_%s_10um.out' %fileIdx
aThread = G4Run(fileName)
aThread.parse()
aRunEventList.extend(aThread.eventList) 
# pickle.dump(aRunEventList, open(fileDir+'/event_list.p', "wb"))

### load psf
psfDir = '/global/scratch/ygzhang/PyBeamDiag/data/psf.mat'
psfTable = []
data_tmp = loadmat(psfDir) 
for i, p in enumerate(data_tmp['psf']):
  psfTable.append(p[30:91, 30:91]) # using 60 x 60 psf

### convert to event object
eventObjList = []
imageTotal = np.zeros((3500, 3500))
numDet = 0
rowHit = []
colHit = []
zHit = []
for event in aRunEventList:
  anEvent = G4Event(event)
  anEvent.detect()
  if anEvent.detected == True:
    numDet += 1
    anEvent.computeHitPosition()
    rowHit.append(anEvent.yHit)
    colHit.append(anEvent.xHit)
    zHit.append(anEvent.zHit)
    anEvent.generateImage(psfTable, pixelPlane='top')
    imageTotal += anEvent.image
  eventObjList.append(anEvent)

colEdges = np.arange(0, 36750.01, 10.5)
rowEdges = colEdges
hitImage, a, b = np.histogram2d(rowHit, colHit, [rowEdges, colEdges])

resultFile = open(fileDir+'/counts.out', 'w')
resultFile.write('%d electrons are produced \n' % len(eventObjList))
resultFile.write('%d electrons are detected on CCD' % numDet)
pickle.dump(hitImage, open(fileDir+'/hitImage_%s.p' % fileIdx, 'wb'))
pickle.dump(imageTotal, open(fileDir+'/trackImage_%s.p' % fileIdx, 'wb'))


'''
save_fits(imageTotal, fileDir+'/10um.fits')
save_fits(hitImage, fileDir+'/hit.fits')
result = {}
result['numDet'] = numDet
result['rowHit'] = rowHit
result['colHit'] = colHit
result['zHit'] = zHit
result['eventObjList'] = eventObjList
result['imageTotal'] = imageTotal
result['hitImage'] = hitImage
pickle.dump(result, open(fileDir+'/process_10um.p', 'wb'))
'''
