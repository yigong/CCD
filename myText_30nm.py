from glob import glob
from os import chdir
G4OutDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/30nm'
chdir(G4OutDir)
fileList = glob('./G4W_*.out')
fId = 0
for i, file in enumerate(fileList):
    f_G4 = open(file, 'r')
    print 'reading ', file
    lines = f_G4.readlines()
    fSplit_name = './split/start_%d.G4track' %(i)
    print 'writing', fSplit_name
    fSplit = open(fSplit_name, 'w')
    for line in lines:
        if line[0] == '*':
            fSplit.close()
            fSplit_name= './split/%d.G4track' %(fId)
            print 'writing', fSplit_name
            fId += 1
            fSplit = open(fSplit_name, 'w')
        elif line == 'Thread-local run terminated.\n':
            fSplit.close()
            fSplit_name = './split/end_%d.G4track' %(i)
            print 'writing', fSplit_name            
            fSplit = open(fSplit_name, 'w')
        else:
            fSplit.write(line)

    fSplit.close()
    f_G4.close()
