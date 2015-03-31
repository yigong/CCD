from CCD.working.Feb24_GammaReconAtExit import EnergyCosPhiParser
from copy import deepcopy

class ElectronEnergyCut:
    ''' A class to apply electron energy cut.
    '''
    

    
    def __init__(self, fileName=None, scatterWidth=10):
        ''' Initialize variables
        '''
        self.fileName = fileName
        self.outputObj = EnergyCosPhiParser(fileName=fileName, scatterWidth=scatterWidth)
        self.outputObj.parse()
        self.outputObj.energyRecon()
        
        cutVars = ('cosPhi', 'energyExit', 'gammaEnergy', 'lineNumber')
        self.outputObj.addCutVar(cutVars)
        self.eCuts = range(0, 1600, 100)
        self.ratioInOverTotal = []
        self.ratioIn = []
        
    def defineEnergyResolution(self, energyResolution=17.):
        self.keep = (self.outputObj.gammaEnergy > (1732 - energyResolution/2)) & \
                    (self.outputObj.gammaEnergy < (1732 + energyResolution/2))
    
    def separateInOut(self):
        self.cutObj_in = deepcopy(self.outputObj)
        self.cutObj_in.cutData(self.keep)
        self.cutObj_out = deepcopy(self.outputObj)
        self.cutObj_out.cutData(~self.keep)
        self.totalEvents = len(self.outputObj.gammaEnergy) * 1.
        
    def cutData(self):
        for eCut in self.eCuts:
            keepIn = self.cutObj_in.energyExit > eCut
            self.cutObj_in.cutData(keepIn)
            print 'keepIn has shape', keepIn.shape
            print 'Nitro_in has shape', self.cutObj_in.gammaEnergy.shape
            keepOut = self.cutObj_out.energyExit > eCut
            self.cutObj_out.cutData(keepOut)
            print 'keepOut has shape', keepOut.shape
            print 'Nitro_out has shape', self.cutObj_out.gammaEnergy.shape
            countsIn_tmp = len(self.cutObj_in.gammaEnergy) * 1.
            self.ratioInOverTotal.append(countsIn_tmp/self.totalEvents)
            countsOut_tmp = len(self.cutObj_out.gammaEnergy) * 1.
            self.ratioIn.append(countsIn_tmp/(countsIn_tmp + countsOut_tmp))
            
    def plot(self, ax, color, styles=['-', '--'], labels=['', '']):
        ax.plot(self.eCuts, self.ratioIn, c=color, linestyle=styles[0], label=labels[0], linewidth=2)
        ax.plot(self.eCuts, self.ratioInOverTotal, c=color, linestyle=styles[1], label=labels[1], linewidth=2)
        
            
        
        