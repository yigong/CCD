import re
from CCD.physics import compton_electron_array as computeElectronEnergy
from mpldatacursor import datacursor
from scipy.optimize import curve_fit
from CCD.gaussian_fit import gaussian_linear_background

mc2 = 511.
BINS = np.linspace(1700, 1800, 100)
binsMid = BINS[:-1] + 0.5*diff(BINS)
# binsFit = np.linspace(1725, 1740, 400)
binsFit = binsMid
figE, axE = plt.subplots()
figEfit, axEfit = plt.subplots()
counts = []
fit_counts = []
Gamma = []
dotOutFiles = ['nitro_10um.out', 'nitro_30um.out', 'nitro_50um.out', 'nitro_70um.out', 'nitro_90um.out', 'nitro_100um.out']
widthStr = ['10 um', '30 um', '50 um', '70 um', '90 um', '100 um']
histShift = [0, 0.01, 0.02, 0.03, 0.04, 0.05]

for fileIdx, fileName in enumerate(dotOutFiles):
    cosPhi = []
    energyExit = []
    lineNumber = []
    process = []
    comptFlag = False
    fileDir = "/Users/Yigong/GEANT4/work/BeamDiag_simple/build/out_GammaReconAtExit/"\
              + fileName 
    fileObj = open(fileDir, "rw+")
    for lineNO, line in enumerate(fileObj):
        if line.find('    0     ')!= -1:
            process_tmp = re.findall(r'\S{1,20}\n', line)
            processName = process_tmp[0][:-1] 
            if processName == 'compt':
                process.append(processName)
                comptFlag = True
        
        
        
        if comptFlag==True and line.find(widthStr[fileIdx])!=-1 and line.find('physiScatter')!=-1 and line.find('Transportation')!=-1:
            energyStr_tmp = re.findall(r'\S{1,15} keV', line)
            energy_tmp = energyStr_tmp[0].split(' ')
            energyExit.append(float(energy_tmp[0]))
            
            dirStr_tmp = re.findall(r'\(.+\)', line)
            dir_tmp = dirStr_tmp[0].split(',')
            cosPhi.append(float(dir_tmp[1]))
            
            lineNumber.append(lineNO)
            comptFlag = False
            #print line
        if line.find('Transportation') != -1:
            comptFlag = False
            
    # for a in zip(energyExit, cosPhi, lineNumber, process):
    #     print a 
    energyExit = np.array(energyExit)
    cosPhi = np.array(cosPhi)
    
    sigmaC = np.sqrt(energyExit*(energyExit+2*mc2))
    gammaEnergy = energyExit*mc2/(cosPhi*sigmaC-energyExit)
    counts_tmp, bins = np.histogram(gammaEnergy, bins=BINS)
    axE.plot(binsMid+histShift[fileIdx], counts_tmp, drawstyle='steps-mid')
    initGuess = [600, 1732, 10., 2.]
    parameters, pCovariance = curve_fit(gaussian_linear_background, binsMid, counts_tmp, initGuess)
    #print parameters
    centroid = parameters[1]
    centroidError = np.sqrt(pCovariance[1,1])
    Std = parameters[-2]
    StdError = np.sqrt(pCovariance[-2,-2])
    FWHM_tmp = Std * 2.35
    peakLow = centroid - 2*Std
    peakHigh = centroid + 2*Std
    label_tmp = '{: <16s}'.format(widthStr[fileIdx]) + '  FWHM = {0:.2f} keV'.format(FWHM_tmp)
    axE.plot(binsMid+histShift[fileIdx], counts_tmp, drawstyle='steps-mid', label=label_tmp)
    counts_fitting = gaussian_linear_background(binsFit, *parameters)
    axEfit.plot(binsFit, counts_fitting, lw=2, label=label_tmp)
    numberInPeak = np.sum(np.logical_and(gammaEnergy>peakLow, gammaEnergy<peakHigh))
    numberElectron = len(gammaEnergy)
    print "Width of {:s}".format(widthStr[fileIdx])
    print "  Centroid: {:.2f} keV  |  FWHM: {:.2f} keV".format(centroid, FWHM_tmp)
    print "  Peak to Total: {:d} out of {:d}".format(numberInPeak, numberElectron)
    
    counts.append(counts_tmp)
    fit_counts.append(counts_fitting)
    Gamma.append(gammaEnergy)
    fileObj.close()
axE.legend()
axEfit.legend()
plt.show()
