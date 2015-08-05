import re
import numpy as np
from scipy.io import loadmat, savemat
from PyBeamDiag.diffuse import XYZdE2track
from astropy.io import fits
from numpy import arctan, rad2deg

def parse(file, psfTable, outFolder, pixelPlane):
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
    xDir = -1
    yDir = -1

    for line in lines:
        if line[0] == '*' and x_tmp:
            print 'see * in ',file 
        elif 'physiTracker' not in line:
            continue
        elif 'physiTracker' in line:
            lineSplit = re.split(r'\s*[(), \s]\s*', line)
            if lineSplit[0] == '':
                lineSplit = lineSplit[1:]
            if not x_tmp:
                xDir = float(lineSplit[7])
                yDir = float(lineSplit[8])
                xInit = float(lineSplit[0]) 
                yInit = float(lineSplit[1])
                eInit = float(lineSplit[5])
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
            alpha_true = rad2deg(arctan(xDir/yDir))
            if alpha_true < 0:
                alpha_true += 360.   
            del x_tmp, y_tmp, z_tmp, dE_tmp
            track, rowMin, colMin = XYZdE2track(x, y, z, dE, psfTable, pixelPlane)
            h = fits.Header()
            h['rowMin'] = rowMin
            h['colMin'] = colMin
            h['alphaT'] = alpha_true
            h['xInit'] = xInit
            h['yInit'] = yInit 
            h['eInit'] = eInit
            fits.writeto('../%s/%s.fits' % (outFolder, fileIdx), track, h, clobber=True)
        f.close()   

def parse_E_phi(g4Out, outFolder):
    '''
    parse electron energy and angle phi
    '''
    f = open(g4Out, 'r')
    lines = f.readlines()
    E_create = []
    cosPhi_create = []
    E_scatterOut = []
    cosPhi_scatterOut = []
    for line in lines:
        if ('3500' in line) and ('compt' in line) and ('physiScatter' in line) :
            lineSplit = re.split(r'\s*[(), \s]\s*', line)
            if lineSplit[0] == '':
                lineSplit = lineSplit[1:]
            E_create.append(lineSplit[5])
            cosPhi_create.append(lineSplit[8])
            
        if ('10' in line) and ('physiScatter' in line) and ('Transportation' in line): 
            lineSplit = re.split(r'\s*[(), \s]\s*', line)
            if lineSplit[0] == '':
                lineSplit = lineSplit[1:]
            E_scatterOut.append(lineSplit[5])
            cosPhi_scatterOut.append(lineSplit[8])

    data_dict = dict()
    data_dict['E_create'] = E_create
    data_dict['cosPhi_create'] = cosPhi_create
    data_dict['E_scatterOut'] = E_scatterOut
    data_dict['cosPhi_scatterOut'] = cosPhi_scatterOut
    savemat(@.mat file name@, data_dict)

