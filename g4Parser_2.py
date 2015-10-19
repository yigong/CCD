import re
import numpy as np
from scipy.io import loadmat, savemat
from PyBeamDiag.diffuse import XYZdE2track
from astropy.io import fits
from numpy import arcsin, arctan2, rad2deg
import os

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
    fileIdx = re.split('\.', file)[0]
    x = []
    y = []
    z = []
    dE = []
    xDir = -1
    yDir = -1

    if os.path.isfile('%s/%s.fits' %(outFolder, fileIdx)):
        return

    for line in f:
        if 'physiTracker' not in line:
            continue
        elif 'physiTracker' in line:
            lineSplit = re.split(r'\s*[(), \s]\s*', line)
            if lineSplit[0] == '':
                lineSplit = lineSplit[1:]
            if not x:
                xDir = float(lineSplit[7])
                yDir = float(lineSplit[8])
                zDir = float(lineSplit[9])
                xInit = float(lineSplit[0]) 
                yInit = float(lineSplit[1])
                zInit = float(lineSplit[2])
                eInit = float(lineSplit[5])
                alphaT = rad2deg(arctan2(-yDir, -xDir))
                betaT = rad2deg(arcsin(-zDir))
            x.append(float(lineSplit[0]))
            y.append(float(lineSplit[1]))
            z.append(float(lineSplit[2]))
            dE.append(float(lineSplit[6]))
    else:
        # executed after 'for' terminates normally
        # when 'for' terminates normally, it reaches the end of file
        if x:
            x = np.array(x)
            y = np.array(y)
            z = np.array(z)
            dE = np.array(dE)
            track, rowEdge0, colEdge0= XYZdE2track(x, y, z, dE, psfTable, pixelPlane)
            h = fits.Header()
            h['row0'] = rowEdge0
            h['col0'] = colEdge0
            h['alphaT'] = alphaT
            h['betaT'] = betaT
            h['xDir'] = xDir
            h['yDir'] = yDir
            h['zDir'] = zDir
            h['xInit'] = xInit
            h['yInit'] = yInit 
            h['zInit'] = zInit
            h['eInit'] = eInit
            prihdu = fits.PrimaryHDU(track, header=h)
            row_pix = (y - rowEdge0)/10.5
            col_pix = (x - colEdge0)/10.5
            row = fits.Column(name='row', format='F', array=row_pix) 
            col = fits.Column(name='col', format='F', array=col_pix)
            tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs([row,col]))
            hdulist = fits.HDUList([prihdu, tbhdu])
            hdulist.writeto('%s/%s.fits' % (outFolder, fileIdx))
        f.close()   

def parse_E_phi(g4Out, outFolder, thk):
    '''
    parse electron energy and angle phi
    '''
    thk = ' ' + thk + ' '
    f = open(g4Out, 'r')
    fNameSplit = re.split(r'[_\.]', g4Out)
    matrl = fNameSplit[2]
    fId   = fNameSplit[1]

    lines = f.readlines()
    E_create = []
    cosPhi_create = []
    E_scatterOut = []
    cosPhi_scatterOut = []
    x_scatterOut = []
    y_scatterOut = []
    z_scatterOut = []
    for line in lines:
        if ('3500' in line) and ('compt' in line) and ('physiScatter' in line) :
            lineSplit = re.split(r'\s*[(), \s]\s*', line)
            if lineSplit[0] == '':
                lineSplit = lineSplit[1:]
            E_create.append(float(lineSplit[5]))
            cosPhi_create.append(float(lineSplit[8]))
            
        if (thk in line) and ('physiScatter' in line) and ('Transportation' in line): 
            lineSplit = re.split(r'\s*[(), \s]\s*', line)
            if lineSplit[0] == '':
                lineSplit = lineSplit[1:]
            E_scatterOut.append(float(lineSplit[5]))
            cosPhi_scatterOut.append(float(lineSplit[8]))
            x_scatterOut.append(float(lineSplit[0]))
            y_scatterOut.append(float(lineSplit[1]))
            z_scatterOut.append(float(lineSplit[2]))
    data_dict = dict()
    data_dict['E_create'] = E_create
    data_dict['cosPhi_create'] = cosPhi_create
    data_dict['E_scatterOut'] = E_scatterOut
    data_dict['cosPhi_scatterOut'] = cosPhi_scatterOut
    data_dict['x_scatterOut'] = x_scatterOut
    data_dict['y_scatterOut'] = y_scatterOut
    data_dict['z_scatterOut'] = z_scatterOut

    savemat('%s/%s_%s.mat' %(outFolder, matrl, fId), data_dict)

