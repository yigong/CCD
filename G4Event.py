import numpy as np
from PyBeamDiag.diffuse import diffuseAndPixelize

class G4Event:
    ''' A class of G4Event
            - converting the output of G4Run to multiple G4Event
            - create a CCD image from 
    '''
    
    def __init__(self, eventDict):
        ''' Initialize parameters
        '''
        
        
        
        self.x = []
        self.y = []
        self.z = []
        self.ds = []
        self.vol = []
        self.E = []
        self.dE = []
        self.dirX = []
        self.dirY = []
        self.dirZ = []
        self.proc = []
        for key in eventDict.keys():
            if hasattr(self, key):
                setattr(self, key, np.array(eventDict[key]))
            else:
                print 'G4Event class doesn\'t have attribute', key
    
    def detect(self):
        ''' set new attributes of 'detected' and 'contained'
        '''
        
        self.inCCDFlag = (self.vol == 'physiTracker')
        self.eTrackX = self.x[self.inCCDFlag]
        self.eTrackY = self.y[self.inCCDFlag]
        self.eTrackZ = self.z[self.inCCDFlag]
        self.eTrackE = self.E[self.inCCDFlag]
        self.eTrackdE = self.dE[self.inCCDFlag][1:]
        self.eTrackds = self.ds[self.inCCDFlag][1:]
        if any(self.inCCDFlag):
            self.detected = True
            if self.inCCDFlag[-1] == False:
                self.contained = False
            else:
                self.contained = False
        else:
            self.detected = False     
            self.contained = False                
        
    def computeHitPosition(self):
        ''' get the hitting position relative to lower left corner in um
        '''
        if self.detected == True:
          self.xHit = self.eTrackX[0] + 37000./2
          self.yHit = self.eTrackY[0] - 2000.
          self.zHit = self.eTrackZ[0]

    def generateImage(self, psfTable, pixelPlane):
        ''' generate a CCD image from a track
        '''
        
        if self.detected == True:
            # convert to CCD coordinate system
            Xtmp = self.eTrackX + 37000./2
            Ytmp = self.eTrackY - 2000.
            Ztmp = self.eTrackZ  
            Etmp = self.eTrackE
            
            # compute the center of each step and energy deposition                
            xCCD = 0.5*(Xtmp[:-1] + Xtmp[1:])
            yCCD = 0.5*(Ytmp[:-1] + Ytmp[1:])
            zCCD = 0.5*(Ztmp[:-1] + Ztmp[1:])
            self.image = diffuseAndPixelize(xCCD, yCCD, zCCD, self.eTrackdE, \
                                            psfTable, pixelPlane)
        else:
            self.image = np.zeros((3500, 3500))
            
            
            
            
            
            
            
            
            
        
