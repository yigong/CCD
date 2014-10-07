from CCD.working.G4Macro import G4Macro
from CCD.working.slurmGenerator import slurmGenerator
from CCD.working.G4Bash import G4Bash

for i in xrange(16):
    g4 = G4Macro('eSample_%s.mac'%(int(i)))
    g4.writeAll(p.energySample[i], p.betaSample[i])
    slm = slurmGenerator(i)
    slm.writeAll()
    
# g4.writePhysics()
# g4.writeSetCuts()
# g4.writeDetectorConstruct()
# g4.writeParticle(energy=1000)
# g4.writeDirection()
# g4.writeBeamOn()

Macros = ['eSample_%s.mac'%(int(i)) for i in range(16)]
bash = G4Bash('eSample.sh')
bash.writeAll(Macros)
