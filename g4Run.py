import re
from IPython.core.debugger import Tracer
import numpy as np
from PyBeamDiag.diffuse import diffuseAndPixelize

class G4Run:
    ''' A class to parse GEANT4 output.
    '''
    def __init__(self, fName):
        self.fileName = fName
        self.endFlag = False 
        f = open(fName, 'r')
        # read the file until the 1st event
        for line in f:
            if line[0] == "*":
                break
        self.f = f
        self.image = np.zeros((3500, 3500))
    
    def parse1Event(self):
        ''' convert an event history to a track on CCD
        ''' 
        x_tmp = []
        y_tmp = []
        z_tmp = []
        dE_tmp = []
        for line in f:
            if line[0] == '*':
                break 
            if 'physiTracker' not in line:
                continue
            if 'physiTracker' in line:
                lineSplit = re.split(r'\s*[(),\s]\s*', line)
                if lineSplit[0] == '':
                    lineSplit = lineSplit[1:]  # remove empty string
                x_tmp.append(float(lineSplit[0]))
                y_tmp.append(float(lineSplit[1]))
                z_tmp.append(float(lineSplit[2]))
                x = np.array(x_tmp) + 37000./2
                y = np.array(y_tmp) - 2000.
                z = np.array(z_tmp)
                dE_tmp.append(float(lineSplit[6]))
                dE = np.array(dE_tmp)
        image_1Track = diffuseAndPixelize(x, y, z, dE, psfTable, pixelPlane='top') 


   #  def __init__(self, fileName):
   #      self.fileName = fileName
   #      self.eventList = []
   #      
   #  def parse(self):
   #      ''' Parse each GEANT4 event to a dict, then append the dict to trackList.
   #      '''
   #      fileObj = open(self.fileName, 'r+')
   #      startFlag = False
   #      anEvent = {}
   #      anEvent['x'] = []
   #      anEvent['y'] = []
   #      anEvent['z'] = []
   #      anEvent['ds'] = []
   #      anEvent['vol'] = []
   #      anEvent['E'] = []
   #      anEvent['dE'] = []
   #      anEvent['dirX'] = []
   #      anEvent['dirY'] = []
   #      anEvent['dirZ'] = []
   #      anEvent['proc'] = []
   #      for line in fileObj:
   #          if line == 'Thread-local run terminated.\n':
   #              self.eventList.pop(0)
   #              startFlag = False

   #          elif startFlag == True and line[0]!='*':
   #              lineSplit = re.split(r'\s*[(),\s]\s*', line)
   #              if lineSplit[0] == '':          
   #                  lineSplit = lineSplit[1:]   # remove the empty string
   #              anEvent['x'].append(float(lineSplit[0]))
   #              anEvent['y'].append(float(lineSplit[1]))
   #              anEvent['z'].append(float(lineSplit[2]))
   #              anEvent['ds'].append(float(lineSplit[3]))
   #              anEvent['vol'].append(lineSplit[4])
   #              anEvent['E'].append(float(lineSplit[5]))
   #              anEvent['dE'].append(float(lineSplit[6]))
   #              anEvent['dirX'].append(float(lineSplit[7]))
   #              anEvent['dirY'].append(float(lineSplit[8]))
   #              anEvent['dirZ'].append(float(lineSplit[9]))
   #              anEvent['proc'].append(lineSplit[10])
   #              
   #          elif startFlag == True and line[0]=='*':
   #              self.eventList.append(anEvent)
   #              anEvent = {}
   #              anEvent['x'] = []
   #              anEvent['y'] = []
   #              anEvent['z'] = []
   #              anEvent['ds'] = []
   #              anEvent['vol'] = []
   #              anEvent['E'] = []
   #              anEvent['dE'] = []
   #              anEvent['dirX'] = []
   #              anEvent['dirY'] = []
   #              anEvent['dirZ'] = []
   #              anEvent['proc'] = []

   #          elif line[0:9] == '### Run 0':
   #              startFlag = True
            

                            
