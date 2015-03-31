from CCD.physics import compton_electron_array as computeElectronEnergy
import re
import numpy as np

class DopplerBroadening:

    mc2 = 511.
    def __init__(self, fileName=None):
        ''' Initialize variables
        '''
        self.cosPhi = []
        self.energyCreation = []
        self.fileName = fileName
       
    def parse(self):
        comptFlag = False
        fileObj = open(self.fileName, "rw+")
        for lineNO, line in enumerate(fileObj):
            if line.find('    0     ')!= -1 and line.find('compt')!=-1 and line.find('3500 um')!=-1 :
                energyStr_tmp = re.findall(r'\S{1,15} keV', line)
                energy_tmp = energyStr_tmp[0].split(' ')
                self.energyCreation.append(float(energy_tmp[0]))
                dirStr_tmp = re.findall(r'\(.+\)', line)
                dir_tmp= dirStr_tmp[0].split(",")
                self.cosPhi.append(float(dir_tmp[1]))
    
        
                
        self.energyCreation = np.array(self.energyCreation)
        self.cosPhi = np.array(self.cosPhi)
        
    def sort(self):
        self.phiRadSorted_tmp, self.G4energySorted_tmp = zip(*sorted(zip(np.arccos(self.cosPhi), self.energyCreation)))
        self.phiRadSorted = np.array(self.phiRadSorted_tmp)
        self.energySorted = np.array(self.G4energySorted_tmp)
        self.phiDegSorted = self.phiRadSorted/np.pi*180.
        self.energyComputed = computeElectronEnergy(203., self.phiRadSorted)
    
    def hist(self, ax, BINS):
        self.energyDiff = self.energySorted - self.energyComputed
        #self.energyDiff = energyDiff
        ax.hist(self.energyDiff, BINS)
        
        