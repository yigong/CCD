import os
import glob
import re
import numpy as np
from scipy.io import loadmat
from PyBeamDiag.io import save_fits
from PyBeamDiag.diffuse import XYZdE2track

# cd to G4out dir
dataDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/10um'
os.chdir(dataDir)
fileList = glob.glob('./G4W_*_10um.out')
# load PSF
psfDir = '/global/scratch/ygzhang/PyBeamDiag/data/psf.mat'
psf_tmp = loadmat(psfDir)
psfTable = np.array(psf_tmp['psf'])
# initialization
nTrk = 0 # count number of tracks
for file in fileList[:2]:
    print 'Now parsing file %s' % file
    endOfFile = False
    f = open(file, 'r')
    while not endOfFile:
        x_tmp = []
        y_tmp = []
        z_tmp = []
        dE_tmp = []
        for line in f:
            if line[0] == '*' and x_tmp:
                break
            elif line == '':
                endOfFile = True             
                print 'Reaches EOF'
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
        x = np.array(x_tmp) + 37000./2
        y = np.array(y_tmp) - 2000.
        z = np.array(z_tmp)
        dE = np.array(dE_tmp)
        track, rowMin, colMin = XYZdE2track(x, y, z, dE, psfTable, \
                pixelPlane='top')
        del x_tmp, y_tmp, z_tmp, dE_tmp, x, y, z, dE
        nTrk += 1
        print 'Now parsing track #%d' % (nTrk)
        print '(rowMin, colMin) = %d, %d' % (rowMin, colMin)
        save_fits(track, './%d.fits' % nTrk, False)
        
