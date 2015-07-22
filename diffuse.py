import numpy as np

#@profile
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
    rowEdges = np.arange(rowMin*10.5, rowMax*10.5 + 0.01, 10.5)
    colEdges = np.arange(colMin*10.5, colMax*10.5 + 0.01, 10.5)
    for i, (x_1step, y_1step, z_1step,dE_1step) in enumerate(zip(xArray, yArray, zArray, dEArray)):
        rows = y_1step + rowPSF1d
        cols = x_1step + colPSF1d
        dE = (dE_1step * computePSF(z_1step, psfTable, pixelPlane))
        track_1step, a, b = np.histogram2d(rows, cols, [rowEdges, colEdges], weights=dE)
        track += track_1step
    del track_1step
    del rows, cols, dE
    return track, rowMin, colMin
    
def computePSF(z, psfTable, pixelPlane):
    ''' compute the point spread function given the z value '''
    if pixelPlane == 'bottom':
        # print 'back-plane irradiation'
        zTable = np.arange(0, 650.5, 0.5)
        zAdjusted = min(zTable, key=lambda x:abs(x-z))        
        return psfTable[int(2*zAdjusted)].flatten()
    
    elif pixelPlane == 'top':
        # print 'pixel-plane irradiation'
        zTable = np.arange(0, 650.5, 0.5)
        zAdjusted = 650 - min(zTable, key=lambda x:abs(x-z))    
        return psfTable[int(2*zAdjusted)].flatten()
    else:
        print 'pixelPlane is neither bottom nor top'
