import os
import glob
import re
import inspect
import numpy as np
from scipy.io import loadmat
from PyBeamDiag.diffuse import XYZdE2track
from astropy.io import fits

# cd to G4out dir
dataDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/3um'
os.chdir(dataDir)
fileList = glob.glob('./G4W_*_3um.out')
# load PSF
psfDir = '/global/scratch/ygzhang/PyBeamDiag/data/psf.mat'
psf_tmp = loadmat(psfDir)
psfTable = np.array(psf_tmp['psf'])
# initialization
lineNum = 0
nTrk = 0 # count number of tracks
for file in fileList:
    print 'Now parsing file %s' % file
    endOfFile = False
    f = open(file, 'r')
    while True:
    # the while loop will be broken when EOF is reached
        try: 
            x_tmp = []
            y_tmp = []
            z_tmp = []
            dE_tmp = []
            for line in f:
                lineNum += 1
                if line[0] == '*' and x_tmp:
                    break
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
                # then break the while-loop
                break

            print 'Now we are at line #%s' % lineNum
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
            h = fits.Header()
            h['rowMin'] = rowMin
            h['colMin'] = colMin
            fits.writeto('./%d.fits' % nTrk, track, h, clobber=True)
        except ValueError:
            print 'ValueError'
            continue
