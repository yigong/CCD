import re

fileObj = open("/Users/Yigong/GEANT4/work/BeamDiag_simple/build/PbCrossSection.out", "rw")
process = []
energy = []
lineNumber = []
gammaFlag = False
for lineNO, line in enumerate(fileObj):
    if line.find('Particle = gamma,   Track ID = 1,   Parent ID = 0')!= -1:
        gammaFlag = True  
    if gammaFlag==True and line.find('    1     ')!= -1:
        process_tmp = re.findall(r'\S{1,20}\n', line)
        process.append(process_tmp[0][:-1])
        energyStr_tmp = re.findall(r'\S{1,15} keV', line)
        energy_tmp = energyStr_tmp[0].split(' ')
        energy.append(float(energy_tmp[0]))
        lineNumber.append(lineNO)
        gammaFlag = False
        
energy = np.array(energy)
process = np.array(process)
lineNumber = np.array(lineNumber)

# cross section data from NIST [cm^2/g]
xSection_Comp = 3.768e-2
xSection_Rayl = 1.014e-3
xSection_PE   = 6.440e-3
xSection_Pair   = 3.490e-3

N_photons = 1e5
areaDensity = 1.135 # in [g/cm^2]


N_Comp_G4 = np.sum(process=='compt')
N_Comp_NIST = N_photons * areaDensity * xSection_Comp
N_Pair_G4 = np.sum(process=='conv')
N_Pair_NIST = N_photons * areaDensity * xSection_Pair
N_PE_G4 = np.sum(process=='phot')
N_PE_NIST = N_photons * areaDensity * xSection_PE
N_Rayl_G4 = np.sum(process=='Rayl')
N_Rayl_NIST = N_photons * areaDensity * xSection_Rayl
print '%s(%s) Compton Scattering' %(N_Comp_G4, int(N_Comp_NIST))
print '%s(%s) Pair Production' %(N_Pair_G4, int(N_Pair_NIST))
print '%s(%s) Photoelectron Absorption' %(N_PE_G4, int(N_PE_NIST))
print '%s(%s) Rayleigh Scattering' %(N_Rayl_G4, int(N_Rayl_NIST))

        