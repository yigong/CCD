import numpy as np

class G4Macro():
    ''' A class to write Macro for GEANT4
    
        An example usage is as follows:
        
        class G4Macro()
        
        then calling:
        G4Macro.
    '''
    
    #
    #
    
    def __init__(self, fileName):
        self.fileObj = open(fileName, 'w')
        
    
    def writeVerbosity(self, controlVerbose=0, runVerbose=0, trackingVerbose=1):
        ''' write the verbosity section '''
        
        lines = ['/control/verbose %s'%(controlVerbose)]
        lines.append('/run/verbose %s'%(runVerbose))
        lines.append('/tracking/verbose %s'%(trackingVerbose))
        lines.append('#\n')
        self.fileObj.write('\n'.join(lines))
        
    def writePhysics(self, physicsList='empenelope'):
        ''' write the PhysicsList section '''
        
        lines = ['/testem/phys/addPhysics %s'%(physicsList)]
        lines.append('#\n')
        self.fileObj.write('\n'.join(lines))
        
    def writeSetCuts(self, cuts='50 nm'):
        ''' write the setCuts section '''
        
        lines = ['/testem/phys/setCuts %s'%(cuts)]
        lines.append('#\n')
        self.fileObj.write('\n'.join(lines))
        
    def writeDetectorConstruct(self, material='Silicon'):
        ''' write the detector construction section ''' 
        
        lines = ['/testem/det/setMat %s'%(material)]
        lines.append('/testem/det/setLbin 1 .00693875')
        lines.append('/testem/det/setRbin 10 .005')
        lines.append('/run/initialize')
        lines.append('/testem/stepMax 100 nm')
        lines.append('#\n')
        self.fileObj.write('\n'.join(lines))
        
    def writeParticle(self, energy, particle='e-'):
        ''' write the particle generator section '''
        
        lines = ['/gun/particle %s'%(particle)]
        lines.append('/gun/energy %.3f keV'%(energy))
        lines.append('#\n')
        self.fileObj.write('\n'.join(lines))
        
    def writeDirection(self, position=(0, 0, -325), direction=(0, 1, 1)):
        ''' write the particle direction section '''
        
        lines = ['/gun/position %s %s %s um'%(position)]
        lines.append('/gun/direction %s %.3f %s um'%(direction))
        lines.append('#\n')
        self.fileObj.write('\n'.join(lines))
        
    def writeBeamOn(self, particleNum=5000):
        ''' write the beamOn line '''
        
        lines = ['/run/beamOn %s'%(particleNum)]
        lines.append('#\n')
        self.fileObj.write('\n'.join(lines))
        self.fileObj.close()
        
    def writeAll(self, energy, beta):
        self.writeVerbosity()
        self.writePhysics()
        self.writeSetCuts()
        self.writeDetectorConstruct()
        self.writeParticle(energy=energy)
        
        yDirection = 1/np.tan(beta)
        direction = (0, yDirection, 1)
        
        self.writeDirection(direction=direction)
        self.writeBeamOn()
        
        
        
        
            
        
    
         
        
        