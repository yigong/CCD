import numpy as np

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
        self.vol = []
        self.E = []
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
        if any(self.inCCDFlag):
            self.detected = True
        if self.inCCDFlag[-1] == False:
            self.contained = False                     
        
    def generateImage(self):
        ''' generate a CCD image from a track
        '''
        
        if self.eTrackX.size != 0:
            # convert to CCD coordinate system
            Xtmp = self.eTrackX + 37000./2
            Ytmp = self.eTrackY - 2000.
            Ztmp = self.eTrackZ  
            Etmp = self.eTrackE
            # compute the center of each step and energy deposition
            xCCD = 0.5*(Xtmp[:-1] + Xtmp[1:])
            yCCD = 0.5*(Ytmp[:-1] + Ytmp[1:])
            zCCD = 0.5*(Ztmp[:-1] + Ztmp[1:])
            dE   = Etmp[:-1] - Etmp[1:]
            
            
            
            
            
            
            
            
            
            
        