from glob import glob
from os import chdir
import sys

def main():

    args = sys.argv
    G4OutDir = args[1]
    splitOutDir = args[2]
    chdir(G4OutDir)
    fileList = glob('./G4W_*.out')
    fId = 0
    for i, file in enumerate(fileList):
        f_G4 = open(file, 'r')
        print 'reading ', file
        lines = f_G4.readlines()
        fSplit_name = '%s/start_%d' %(splitOutDir, i)
        print 'writing ', fSplit_name
        fSplit = open(fSplit_name, 'w')
        for line in lines:
            if line[0] == '*':
                fSplit.close()
                fSplit_name= '%s/%d.G4track' %(splitOutDir, fId)
                print 'writing ', fSplit_name
                fId += 1
                fSplit = open(fSplit_name, 'w')
            elif line == 'Thread-local run terminated.\n':
                fSplit.close()
                fSplit_name = '%s/end_%d' %(splitOutDir, i)
                print 'writing ', fSplit_name            
                fSplit = open(fSplit_name, 'w')
            else:
                fSplit.write(line)

        fSplit.close()
        f_G4.close()

if __name__=='__main__':
    main()

