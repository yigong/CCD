import numpy as np

@profile
def XYZdE2track(xArray, yArray, zArray, dEArray, psfTable, pixelPlane):
    ''' compute the diffused energy deposition'''

    # find minimum and maximum locations in x and y
    xMin = min(xArray)
    xMax = max(xArray)
    yMin = min(yArray)
    yMax = max(yArray)
    # find the row number and col number to use in the track
    rowMin = np.ceil(yMin/10.5) - 3
    rowMax = np.ceil(yMax/10.5) + 3
    colMin = np.ceil(xMin/10.5) - 3
    colMax = np.ceil(xMax/10.5) + 3
    nrow = rowMax - rowMin
    ncol = colMax - colMin
    # create a track array: nrow x ncol
    track = np.zeros((nrow, ncol))
    # compute PSF grid
    rowPSF = np.arange(-30, 30.01, 0.5)
    colPSF = rowPSF
    rowPSF2d, colPSF2d = np.meshgrid(rowPSF, colPSF, indexing='ij')
    rowPSF1d = rowPSF2d.flatten()
    colPSF1d = colPSF2d.flatten()
    # store data
    rows = []
    cols = []
    dE = [] 
    for x_1step, y_1step, z_1step,dE_1step in zip(xArray, yArray, zArray, dEArray):
        rows.extend(y_1step + rowPSF1d)
        cols.extend(x_1step + colPSF1d)
        dE.extend(dE_1step * computePSF(z_1step, psfTable, pixelPlane='top'))
    rowEdges = np.arange(rowMin*10.5, rowMax*10.5 + 0.01, 10.5)
    colEdges = np.arange(colMin*10.5, colMax*10.5 + 0.01, 10.5)
    track, a, b = np.histogram2d(rows, cols, [rowEdges, colEdges], weights=dE) 
    del rows, cols, dE
    return track, rowMin, colMin
##     nrow = psfTable[0].shape[0]
##     ncol = nrow
##     rows = np.arange(-(nrow/2), nrow/2., 1)
##     cols = rows
##     rowVectorPSF, colVectorPSF = np.meshgrid(rows, cols, indexing='ij')
##     # debug_here()
##     rowData = []
##     colData = []
##     EData = []
##     
##     for x, y, z, dE in zip(xArray, yArray, zArray, dEArray):
##         dEv = dE * computePSF(z, psfTable, pixelPlane)
##         rowVector = y + rowVectorPSF
##         colVector = x + colVectorPSF
##         
##         rowData.extend(rowVector.flatten())
##         colData.extend(colVector.flatten())
##         EData.extend(dEv.flatten())
##         # debug_here()
##     colEdges = np.arange(0, 36750.01, 10.5)
##     rowEdges = colEdges
##     image, a, b = np.histogram2d(rowData, colData, [colEdges, rowEdges], weights=EData)
##     # the first element is row; the second is col
##     return image
    
    
def computePSF(z, psfTable, pixelPlane):
    ''' compute the point spread function given the z value '''
    if pixelPlane == 'bottom':
        zTable = np.arange(0, 650.5, 0.5)
        zAdjusted = min(zTable, key=lambda x:abs(x-z))        
        return psfTable[int(2*zAdjusted)].flatten()
    
    elif pixelPlane == 'top':
        zTable = np.arange(0, 650.5, 0.5)
        zAdjusted = 650 - min(zTable, key=lambda x:abs(x-z))    
        return psfTable[int(2*zAdjusted)].flatten()
    else:
        print 'pixelPlane is neither bottom nor top'
