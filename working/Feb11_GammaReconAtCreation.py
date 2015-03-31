import re
from CCD.physics import compton_electron_array as computeElectronEnergy
#from mpldatacursor import datacursor
from scipy.optimize import curve_fit
from CCD.gaussian_fit import gaussian_no_background

outputFiles = ['lithium_21mm.out', 'nitro_7mm.out', 'lead_1mm.out']
scatterStr = ['Lithium', 'Nitrocellulose', 'Lead']
scatterClr = ['b', 'r', 'y']
scatterLineS = ['-', '-', '-']
histShift = [0, 0.02, 0.04]
mc2 = 511.
BINS = np.arange(1720, 1745, 0.1)
binsMid = BINS[:-1] + 0.5*diff(BINS)
binsFit = np.arange(1720, 1745, 0.01)

figE, axE = plt.subplots(figsize=(16, 8))

for i, fileName in enumerate(outputFiles):
    fileDir = "/Users/Yigong/GEANT4/work/BeamDiag_simple/build/out_GammaReconAtCreation/"\
              + fileName 
    fileObj = open(fileDir, "rw+")
    cosPhi = []
    G4energy = []
    eFlag = False
    for line in fileObj:
        #print line
        if line.find('Particle = e-,   Track ID = 2,   Parent ID = 1')!= -1:
            eFlag = True
        if line.find('Particle = gamma,   Track ID = 1,   Parent ID = 0')!= -1:
            eFlag = False
        if eFlag == True and (line.find("compt")!= -1):
            energyStr_tmp = re.findall(r'\S{1,15} keV', line)
            energy_tmp = energyStr_tmp[0].split(' ')
            G4energy.append(float(energy_tmp[0]))
            dirStr_tmp = re.findall(r'\(.+\)', line)
            dir_tmp= dirStr_tmp[0].split(",")
            cosPhi.append(float(dir_tmp[1]))
            eFlag = False
        if (line.find("Track ID = 3")!= -1):
            eFlag = False
            
    phiRad = np.arccos(cosPhi)
    phiRadSorted_tmp, G4energySorted_tmp = zip(*sorted(zip(phiRad, G4energy)))
    phiRadSorted = np.array(phiRadSorted_tmp)
    G4energySorted = np.array(G4energySorted_tmp)
    phiDegSorted = phiRadSorted/np.pi*180.
    ComputedElectronEnergy = computeElectronEnergy(1732., phiRadSorted)
#     fig, ax = plt.subplots()
#     ax.plot(phiDegSorted, ComputedElectronEnergy, 'r', label='Compton Kinematics')
#     ax.plot(phiDegSorted, G4energySorted, color='blue', ls='', alpha=0.3, marker='o', markersize=4, markeredgewidth=0, label='GEANT4')
#     ax.legend()
#     ax.set_xlabel('Electron Scattering Angle ($^{\circ}$)')
#     ax.set_ylabel('Electron Energy (keV)')
#     # ax.set_xlim([0, 90])
#     ax.set_title(scatterStr[i])
    nEvent = len(G4energySorted) 
    print "nEvent when using {}:".format(scatterStr[i]), nEvent
    cosPhi = np.cos(phiRadSorted)
    sigmaC = np.sqrt(G4energySorted*(G4energySorted+2*mc2))
    gammaEnergy = G4energySorted*mc2/(cosPhi*sigmaC-G4energySorted)
    counts_tmp, bins = np.histogram(gammaEnergy, bins=BINS)
    counts_norm = counts_tmp
    
    initGuess = [10., 1732., 2.]
    parameters, pCovariance = curve_fit(gaussian_no_background, binsMid, counts_norm, initGuess)
    FWHM_tmp = parameters[-1] * 2.35
    label_tmp = '{: <16s}'.format(scatterStr[i]) + '  FWHM = {0:.2f} keV'.format(FWHM_tmp)
    counts_fitting = gaussian_no_background(binsFit, *parameters)
    axE.plot(binsMid+histShift[i], counts_norm, c=scatterClr[i], lw=1, ls=scatterLineS[i],\
            drawstyle='steps-mid', label=label_tmp)
    axE.plot(binsFit, counts_fitting, c=scatterClr[i], lw=2,)
    print parameters
#axE.grid()
axE.legend(loc=2)
axE.set_ylabel('Counts')
axE.set_xlabel('Gamma-ray Energy (keV)    binWidth = {} keV'.format(0.1))
axE.set_xlim(1720, 1745)
#datacursor()

plt.show()