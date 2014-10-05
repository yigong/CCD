import numpy as np
from CCD.calibrate import self_subtraction


class Spectrum:
    ''' A class to get the spectrum of a data set
    
    '''
    
    def __init__(self, files):
        self.files = files
        
    def correct(self):
        for file in self.files:
            self_subtraction(file, save_flag=True)
            
        
    def segment(self):
                
        