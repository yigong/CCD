import numpy as np

class PositionSampling():
    ''' A class to sample positions on CCD
        
        An example usage is as follows:
        
        class PositionSampling(alphaMap, betaMap, phiMap, energyMap)
        
        then calling:
        PositionSampling.sample()
        We will get samples on alpha, beta and phi
    '''
    
    def __init__(self, alphaMap, betaMap, phiMap, energyMap, numPoints=16):
        self.numPoints1d = np.sqrt(numPoints)
        self.alphaMap = alphaMap
        self.betaMap = betaMap
        self.phiMap = phiMap
        self.energyMap = energyMap
    
    def sample(self):
        ''' Sample on alpha, beta, phi and energy '''
        
        index1D = np.linspace(0, 3505, self.numPoints1d+2)[1:-1]
        xSampleVector, ySampleVector = np.meshgrid(index1D, index1D, indexing='ij')
        xSample = xSampleVector.flatten().astype(int)
        ySample = ySampleVector.flatten().astype(int)
        
        self.xSample = xSample
        self.ySample = ySample
        self.alphaSample = self.alphaMap[xSample, ySample]
        self.betaSample = self.betaMap[xSample, ySample]
        self.phiSample = self.phiMap[xSample, ySample]
        self.energySample = self.energyMap[xSample, ySample]
        

