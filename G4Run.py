import re
from IPython.core.debugger import Tracer

class G4Run:
    ''' A class to parse GEANT4 output.
    '''
    
    def __init__(self, fileName):
        self.fileName = fileName
        self.eventList = []
        
    def parse(self):
        ''' Parse each GEANT4 event to a dict, then append the dict to trackList.
        '''
        fileObj = open(self.fileName, 'r+')
        startFlag = False
        anEvent = {}
        anEvent['x'] = []
        anEvent['y'] = []
        anEvent['z'] = []
        anEvent['vol'] = []
        anEvent['E'] = []
        anEvent['dirX'] = []
        anEvent['dirY'] = []
        anEvent['dirZ'] = []
        anEvent['proc'] = []
        Tracer()
        for line in fileObj:
            if line == 'Run terminated.\n':
                startFlag = False

            elif startFlag == True and line !='\n':
                lineSplit = re.split(r'\s*[(),\s]\s*', line)
                if lineSplit[0] == '':          
                    lineSplit = lineSplit[1:]   # remove the empty string
                anEvent['x'].append(float(lineSplit[0]))
                anEvent['y'].append(float(lineSplit[1]))
                anEvent['z'].append(float(lineSplit[2]))
                anEvent['vol'].append(lineSplit[3])
                anEvent['E'].append(float(lineSplit[4]))
                anEvent['dirX'].append(float(lineSplit[5]))
                anEvent['dirY'].append(float(lineSplit[6]))
                anEvent['dirZ'].append(float(lineSplit[7]))
                anEvent['proc'].append(lineSplit[8])
                
            elif startFlag == True and line=='\n':
                self.eventList.append(anEvent)
                anEvent = {}
                anEvent['x'] = []
                anEvent['y'] = []
                anEvent['z'] = []
                anEvent['vol'] = []
                anEvent['E'] = []
                anEvent['dirX'] = []
                anEvent['dirY'] = []
                anEvent['dirZ'] = []
                anEvent['proc'] = []

            elif line == '### Run 0 starts.\n':
                startFlag = True
            

                            