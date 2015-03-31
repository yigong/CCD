import re
from CCD.physics import compton_electron_array as computeElectronEnergy
#from mpldatacursor import datacursor
from scipy.optimize import curve_fit
from CCD.gaussian_fit import gaussian_no_background
from CCD.gaussian_fit import doubleGaussianFit
import numpy as np

class EnergyCosPhiParser:
    ''' A class to parse electron energy and direction from GEANT4 output.
    '''
    mc2 = 511.
    
    def __init__(self, fileName=None, scatterWidth=10):
        ''' Initialize variables
        '''
        self.cutVars = []
        self.cosPhi = []
        self.cosPhiCreation = []
        self.energyExit = []
        self.energyCreation = []
        self.lineNumber = []
        self.process = []
        self.fileName = fileName
        self.widthStr = '{:.0f} um'.format(scatterWidth)

    def parse(self):
        energyCreation_tmp = []
        cosPhiCreation_tmp = []
        comptFlag = False
        fileObj = open(self.fileName, "rw+")
        for lineNO, line in enumerate(fileObj):
            if line.find('    0     ')!= -1:
                process_tmp = re.findall(r'\S{1,20}\n', line)
                processName = process_tmp[0][:-1]
                if processName == 'compt':
                    self.process.append(processName)
                    comptFlag = True
            
            if line.find('    0     ')!= -1 and line.find('compt')!=-1 and line.find('3500 um')!=-1 :
                energyStr_tmp = re.findall(r'\S{1,15} keV', line)
                energy_tmp = energyStr_tmp[0].split(' ')
                energyCreation_tmp.append(float(energy_tmp[0]))
                            
                dirStr_tmp = re.findall(r'\(.+\)', line)
                dir_tmp = dirStr_tmp[0].split(',')
                cosPhiCreation_tmp.append(float(dir_tmp[1]))
                
                
            if comptFlag==True and line.find(self.widthStr)!=-1\
               and line.find('physiScatter')!=-1 and line.find('Transportation')!=-1:
                energyStr_tmp = re.findall(r'\S{1,15} keV', line)
                energy_tmp = energyStr_tmp[0].split(' ')
                self.energyExit.append(float(energy_tmp[0]))
            
                dirStr_tmp = re.findall(r'\(.+\)', line)
                dir_tmp = dirStr_tmp[0].split(',')
                self.cosPhi.append(float(dir_tmp[1]))
            
                self.energyCreation.append(energyCreation_tmp[-1])
                self.cosPhiCreation.append(cosPhiCreation_tmp[-1])
                self.lineNumber.append(lineNO)
                comptFlag = False
            
            if line.find('Transportation') != -1:
                comptFlag = False
                
        self.energyExit = np.array(self.energyExit)
        self.cosPhi = np.array(self.cosPhi)
        self.lineNumber = np.array(self.lineNumber)
        self.process = np.array(self.process)
        self.cosPhiCreation = np.array(self.cosPhiCreation)
        self.energyCreation = np.array(self.energyCreation)
    def computeEnergyAngleDiff(self):
        self.angleDiff = np.arccos(self.cosPhiCreation) - np.arccos(self.cosPhi)
        self.energyDiff = self.energyCreation - self.energyExit 
    
    def energyRecon(self):
        sigmaC = np.sqrt(self.energyExit*(self.energyExit+2*self.mc2))
        self.gammaEnergy = self.energyExit*self.mc2/(self.cosPhi*sigmaC-self.energyExit)
    
    def addCutVar(self, attr):
        ''' add attributes to be cut'''
        self.cutVars.extend(attr)
    
    def cutData(self, keep):
        print 'Cutting data from {:d} to {:d}'.format(keep.shape[0], keep.sum())
        for var in self.cutVars:
            if hasattr(self, var):
                ar = getattr(self, var)
                if ar.shape[0] == keep.shape[0]:
                    ar = ar[keep]
                    setattr(self, var, ar)
                else:
                    print '{:s} not same shape: {:s}'.format(var, ar.shape)

            else:
                print '{:s} not found'.format(var)
                
    def printData(self):
        print 'Variable name (shape):'
        for var in self.cutVars:
            if hasattr(self, var):
                ar = getattr(self, var)
                print '{:s}, {:s}'.format(var, ar.shape)
            else:
                print '{:s} not found'.format(var)
    
    def hist(self, ax, color='b', bins=np.arange(1600,1800,1), norm=False):
        counts_tmp, bins = np.histogram(self.gammaEnergy, bins=bins)
        self.counts = counts_tmp
        if norm == True:
            counts = counts_tmp.astype(float)/counts_tmp.max()
        else:
            counts = counts_tmp
        self.binsMid = bins[:-1] + 0.5*np.diff(bins)
        self.binWidth = bins[1] - bins[0]
        ax.plot(self.binsMid, counts, label=self.widthStr, drawstyle='steps-mid', color=color)
        ax.legend()
    def GaussianFit(self, ax, color='b', initGuess=[1000, 1732, 10.]):
        parameters, pCovariance = curve_fit(gaussian_no_background, \
                                            self.binsMid, self.counts, initGuess)
        print self.widthStr, parameters
        centroid = parameters[1]
        centrqoidError = np.sqrt(pCovariance[1,1])
        Std = parameters[2]
        StdError = np.sqrt(pCovariance[2,2])
        FWHM_tmp = Std*2.35
        peakLow = centroid - 2*Std
        peakHigh = centroid + 2*Std
        label_tmp = self.widthStr + '   FWHM = {0:.2f} keV'.format(FWHM_tmp)
        binsFit = np.linspace(1600, 1800, 1000)
        counts_fitting = gaussian_no_background(binsFit, *parameters)
        ax.plot(binsFit, counts_fitting, label=label_tmp, color=color)
        ax.legend()
        
        
