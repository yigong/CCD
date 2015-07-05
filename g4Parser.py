import re
import numpy as np
from scipy.io import loadmat
from PyBeamDiag.diffuse import XYZdE2track
from astropy.io import fits

#@profile
def parse(psfTable, file):
    '''
    a parser for G4output.
    It's easier to parallelize if written as a separate function.
    '''
    # initialization
    lineNum = 0
    nTrk = 0 # count number of tracks
    # prepare files 
    print 'parsing # ', file
    f = open(file, 'r')
    lines = f.readlines()
    fileIdx = re.split('\.', file)[0]
    x_tmp = []
    y_tmp = []
    z_tmp = []
    dE_tmp = []
    for line in lines:
        if line[0] == '*' and x_tmp:
            print 'see * in ',file 
        elif 'physiTracker' not in line:
            continue
        elif 'physiTracker' in line:
            lineSplit = re.split(r'\s*[(), \s]\s*', line)
            if lineSplit[0] == '':
                lineSplit = lineSplit[1:]
            x_tmp.append(float(lineSplit[0]))
            y_tmp.append(float(lineSplit[1]))
            z_tmp.append(float(lineSplit[2]))
            dE_tmp.append(float(lineSplit[6]))
    else:
        # executed after 'for' terminates normally
        # when 'for' terminates normally, it reaches the end of file
        if x_tmp:
            x = []
            y = []
            z = []
            dE = []
            x = np.array(x_tmp) + 37000./2
            y = np.array(y_tmp) - 2000.
            z = np.array(z_tmp)
            dE = np.array(dE_tmp)
            del x_tmp, y_tmp, z_tmp, dE_tmp
            track, rowMin, colMin = XYZdE2track(x, y, z, dE, psfTable, \
                    pixelPlane='top')
            h = fits.Header()
            h['rowMin'] = rowMin
            h['colMin'] = colMin
            fits.writeto('../fits/%s.fits' % (fileIdx), track, h, clobber=True)
        f.close()   
