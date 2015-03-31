import numpy as np
from IPython.core.debugger import Tracer
debug_here = Tracer()

def diffuseAndPixelize(xArray, yArray, zArray, dEArray, psfTable, pixelPlane):
    ''' compute the diffused energy deposition'''
    nrow = psfTable[0].shape[0]
    ncol = nrow
    rows = np.arange(-(nrow/2), nrow/2., 1)
    cols = rows
    rowVectorPSF, colVectorPSF = np.meshgrid(rows, cols, indexing='ij')
    # debug_here()
    rowData = []
    colData = []
    EData = []
    
    for x, y, z, dE in zip(xArray, yArray, zArray, dEArray):
        dEv = dE * computePSF(z, psfTable, pixelPlane)
        rowVector = y + rowVectorPSF
        colVector = x + colVectorPSF
        
        rowData.extend(rowVector.flatten())
        colData.extend(colVector.flatten())
        EData.extend(dEv.flatten())
        # debug_here()
    colEdges = np.arange(0, 36750.01, 10.5)
    rowEdges = colEdges
    image, a, b = np.histogram2d(colData, rowData, [colEdges, rowEdges], weights=EData)
    
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
