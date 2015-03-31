import re
from CCD.physics import compton_electron_array as computeElectronEnergy

fileObj = open("/Users/Yigong/GEANT4/work/BeamDiag_simple/build/mylarPene_14um.out", "rw+")
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
fig, ax = plt.subplots()
ax.plot(phiDegSorted, ComputedElectronEnergy, 'r', label='Compton Kinematics')
ax.plot(phiDegSorted, G4energySorted, color='blue', ls='', alpha=0.3, marker='o', markersize=4, markeredgewidth=0, label='GEANT4')
ax.legend()
ax.set_xlabel('Electron Scattering Angle ($^{\circ}$)')
ax.set_ylabel('Electron Energy (keV)')
ax.set_xlim([0, 90])
ax.set_title('Lead', fontsize='xx-large')
diff_Pb = G4energySorted - ComputedElectronEnergy

 
fileObj = open("/Users/Yigong/GEANT4/work/BeamDiag_simple/build/LiPene.out", "rw+")
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
fig, ax = plt.subplots()
ax.plot(phiDegSorted, ComputedElectronEnergy, 'r', label='Compton Kinematics')
ax.plot(phiDegSorted, G4energySorted, color='blue', ls='', alpha=0.3, marker='o', markersize=4, markeredgewidth=0, label='GEANT4')
ax.legend()
ax.set_xlabel('Electron Scattering Angle ($^{\circ}$)')
ax.set_ylabel('Electron Energy (keV)')
ax.set_xlim([0, 90])
ax.set_title('Lithium', fontsize='xx-large')
diff_Li = G4energySorted - ComputedElectronEnergy
