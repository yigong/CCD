import numpy as np

class G4Bash():
    ''' a class to write a bash script for job submission on LRC '''
    
    def __init__(self, fileName):
        
        self.fileObj = open(fileName, 'w')
        
    def writeAll(self, Macros):
        lines = []
        for Macro in Macros:
            lines.append('sbatch %s'%(Macro))
        self.fileObj.write('\n'.join(lines))
        self.fileObj.close()
            
        
        