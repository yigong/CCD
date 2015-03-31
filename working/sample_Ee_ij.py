import numpy as np
from copy import copy
from CCD.physics import compton_electron_array as computeEe
from numpy.random import random_sample
from scipy.optimize import curve_fit

colDiameter = 3.
numPhi = 100000

r_sampled = colDiameter/2 * random_sample(numPhi)**0.5 # r = epsilon**(1/2)
theta_sampled = 2 * np.pi * random_sample(numPhi) 
windowX_sampled = r_sampled * np.cos(theta_sampled)
windowZ_sampled = r_sampled * np.sin(theta_sampled) + colDiameter/2 + 1

rowsTmp = np.linspace(0, 3505, 5).astype(int)
colsTmp = np.linspace(0, 3505, 5).astype(int)
[rows, cols] = np.meshgrid(rowsTmp, colsTmp, indexing='ij')
rows = rows.flatten()
cols = cols.flatten()

Egamma = []
# define gaussian fitting function
def GaussFunc(x, totalArea, meanVal, std):
    return totalArea/((2*np.pi)**0.5 * std) * np.exp(-(x - meanVal)**2/(2 * std**2))

for i, (rowPixel, colPixel) in enumerate(zip(rows, cols)): 
    yEdge = 2.97
    zEdge = 1.0
    pxSize = 10.5/1000 # in mm
    mc2 = 511.
    gammaEnergy = 1732.
    CCDlowerLeft = {"x": -1753.5*pxSize, "y": yEdge+0.5*pxSize, "z": 0} # in mm
    
    xPixel = CCDlowerLeft["x"] + colPixel * pxSize
    yPixel = CCDlowerLeft["y"] + rowPixel * pxSize
    zPixel = 0
    alphaRadius = np.abs(np.arctan2(yPixel, xPixel) - np.pi/2)
    cosAlpha = np.cos(alphaRadius)
    
    windowLower = {"x": 0, "y": 0, "z": zEdge} # in mm
    xWindow = windowX_sampled + windowLower["x"]
    yWindow = np.zeros_like(xWindow)
    zWindow = windowZ_sampled + windowLower["z"]
    cosBeta = np.cos(np.arctan2(zWindow, np.sqrt(xPixel**2 + yPixel**2)))
    
    cosPhi = yPixel/np.sqrt((xPixel-xWindow)**2 + \
                            (yPixel-yWindow)**2 + (zPixel-zWindow)**2)
    PhiRadius = np.arccos(cosPhi)
    Ee = computeEe(gammaEnergy, PhiRadius)
    
    numEgamma = 100000
    counts_Ee, binEdges_Ee = np.histogram(Ee, bins=100, normed=True)
    binWidth_Ee = binEdges_Ee[1] - binEdges_Ee[0] 
    bins_Ee = binEdges_Ee[:-1] + binWidth_Ee/2.
    cdf_Ee = np.cumsum(counts_Ee*binWidth_Ee)
    EeSampled = np.interp(random_sample(numEgamma), cdf_Ee, bins_Ee)
    
    counts_Phi, binEdges_Phi = np.histogram(PhiRadius, bins=100, normed=True)
    binWidth_Phi = binEdges_Phi[1] - binEdges_Phi[0]
    bins_Phi = binEdges_Phi[:-1] + binWidth_Phi/2
    cdf_Phi = np.cumsum(counts_Phi*binWidth_Phi)
    PhiSampled = np.interp(random_sample(numEgamma), cdf_Phi, bins_Phi)
    
    cosPhiSampled = np.cos(PhiSampled)
    sigmaC = (EeSampled * (EeSampled + 2 * mc2))**0.5
    EgammaTmp = EeSampled * mc2 / (cosPhiSampled * sigmaC - EeSampled)
    Egamma.append(EgammaTmp)
    fig, ax = plt.subplots()
    fig.subplots_adjust(bottom=0.13)
    bins = linspace(1000, 2000, 201)
    (counts, binEdges, dummy) = ax.hist(EgammaTmp, bins)
    p0 = [np.sum(counts), np.mean(binEdges), 50]
    fitParameter, fitCovariances = curve_fit(GaussFunc, binEdges[:-1], counts, p0)
    ax.set_title('FWHM = %.2f %% @ i,j = (%d,%d) ' %(2.35*fitParameter[2]/1732*100., rowPixel, colPixel))
    ax.plot(binEdges, GaussFunc(binEdges, *fitParameter), 'r')
    ax.set_xlabel('$E_{\gamma}$ (keV)', fontsize='xx-large')
    ax.set_ylabel('Counts')
    ax.vlines(1732., 0, np.max(counts), color='r', lw=2 )
    
plt.show()
