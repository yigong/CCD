from CCD.G4Run import G4Run
from CCD.G4Event import G4Event
from CCD.io import save_fits
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt

# initialize a run
fileName = '/Users/Yigong/GEANT4/work/BeamDiag_simple/build_1/10um.out'
aRun = G4Run(fileName)
aRun.parse()
aRun.eventList.pop(0)

# load psf
psfTable = []
data = loadmat('psf.mat')
for i, p in enumerate(data['psf']):
    psfTable.append(p[30:91, 30:91])

# convert to event object    
eventObjList = []
imageTotal = np.zeros((3500, 3500))
numDet = 0
for event in aRun.eventList:
    anEvent = G4Event(event)
    anEvent.detect()
    if anEvent.detected == True:
        numDet += 1
        anEvent.generateImage(psfTable, pixelPlane='top')
        imageTotal += anEvent.image
    eventObjList.append(anEvent)

# plot
print '%d Compton electrons are detected.'
save_fits(imageTotal, '10um_image.fits', True)
fig, ax = plt.subplots()
im = ax.imshow(imageTotal, origin='lower')
cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
fig.colorbar(im, cax=cax)
plt.show()
    