from CCD.G4Event import G4Event
# load psf
from scipy.io import loadmat
psfTable = []
data = loadmat('psf.mat')
for i, p in enumerate(data['psf']):
    psfTable.append(p[30:91, 30:91])
    
eventObjList = []
imageTotal = np.zeros((3500, 3500))

for event in aRun.eventList:
    anEvent = G4Event(event)
    anEvent.detect()
    if anEvent.detected == True:
        anEvent.generateImage(psfTable, pixelPlane='top')
        imageTotal += anEvent.image
    eventObjList.append(anEvent)
    
fig, ax = plt.subplots()
im = ax.imshow(imageTotal, origin='lower')
cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
fig.colorbar(im, cax=cax)
plt.show()
    