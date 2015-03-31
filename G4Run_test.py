from CCD.G4Run import G4Run

fileName = '/Users/Yigong/GEANT4/work/BeamDiag_simple/build_Xcode/a.out'
aRun = G4Run(fileName)
aRun.parse()
aRun.eventList.pop(0)