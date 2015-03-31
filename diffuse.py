import numpy as np
from IPython.core.debugger import Tracer
debug_here = Tracer()

def diffuseAndPixelize(xArray, yArray, zArray, dEArray, psfTable, pixelPlane):
    ''' compute the diffused energy deposition'''
    nx = psfTable[0].shape[0]
    ny = nx
    x = np.arange(-(nx/2), nx/2., 1)
    y = x
    xv_tmp, yv_tmp = np.meshgrid(x, y)
    xData = []
    yData = []
    EData = []
    
    for x, y, z, dE in zip(xArray, yArray, zArray, dEArray):
        dEv = dE * computePSF(z, psfTable, pixelPlane)
        xv = x + xv_tmp
        yv = y + yv_tmp
        xData.extend(xv.flatten())
        yData.extend(yv.flatten())
        EData.extend(dEv.flatten())
        debug_here()
    xEdges = np.arange(0, 36750.01, 10.5)
    yEdges = xEdges
    image, a, b = np.histogram2d(xData, yData, [xEdges, yEdges], weights=EData)
    
    return image
    
    
def computePSF(z, psfTable, pixelPlane):
    ''' compute the point spread function given the z value '''
    if pixelPlane == 'bottom':
        zTable = np.arange(0, 650.5, 0.5)
        zAdjusted = min(zTable, key=lambda x:abs(x-z))        
        return psfTable[int(2*zAdjusted)]
    
    elif pixelPlane == 'top':
        zTable = np.arange(0, 650.5, 0.5)
        zAdjusted = 650 - min(zTable, key=lambda x:abs(x-z))    
        return psfTable[int(2*zAdjusted)]
    else:
        print 'pixelPlane is neither bottom nor top'
