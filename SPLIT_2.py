from glob import glob
from os import chdir
import sys
import re
import traceback

try:
    args = sys.argv
    G4Dir = args[1]
    splitDir = args[2]
    chdir(G4Dir)
    fileList = glob('./G4W_*.out')
    fId = 0

    for i, file in enumerate(fileList):
        fG4 = open(file, 'r')

        for line in fG4:
            if re.search('###', line):
                fSplit = []
                break

        for line in fG4:
            if re.search('NewEvent', line):
                if fSplit:              # check if fSplit opened a file
                    fSplit.close()
                fSplit = open('%s/%d.G4track' %(splitDir, fId), 'w')
                fId += 1
            elif re.search('Thread-local run terminated', line):
                fSplit.close()
                break
            else:
                fSplit.write(line)

        fG4.close()
except:
    traceback.print_exc()


