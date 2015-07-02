import re
import numpy as np
from scipy.io import loadmat
from PyBeamDiag.diffuse import XYZdE2track
from astropy.io import fits
import gc

@profile
def parse(file):
    '''
    a parser for G4output.
    It's easier to parallelize if written as a separate function.
    '''
    # initialization
    lineNum = 0
    nTrk = 0 # count number of tracks
    # load PSF
    psfDir = '/global/scratch/ygzhang/PyBeamDiag/data/psf.mat'
    psf_tmp = loadmat(psfDir)
    psfTable = np.array(psf_tmp['psf'])
    del psf_tmp
    # prepare files 
    f = open(file, 'r')
    fileIdx = re.split('_', file)[1]
    fLog = open('./W%s_parser.log' % fileIdx, 'w')
    fLog.write('Now parsing file %s\n' % file)
    while True:
    # the while loop will be broken when EOF is reached
        try: 
            x_tmp = []
            y_tmp = []
            z_tmp = []
            dE_tmp = []
            x = []
            y = []
            z = []
            dE = []
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

            fLog.write('Now we are at line #%s\n' % lineNum)
            x = np.array(x_tmp) + 37000./2
            y = np.array(y_tmp) - 2000.
            z = np.array(z_tmp)
            dE = np.array(dE_tmp)
            track, rowMin, colMin = XYZdE2track(x, y, z, dE, psfTable, \
                    pixelPlane='top')
            nTrk += 1
            fLog.write('Now parsing track #%d\n' % (nTrk))
            fLog.write('(rowMin, colMin) = %d, %d\n' % (rowMin, colMin))
            h = fits.Header()
            h['rowMin'] = rowMin
            h['colMin'] = colMin
            fits.writeto('./fits/W%s_%d.fits' % (fileIdx, nTrk), track, h, clobber=True)
        except ValueError:
            print 'ValueError\n'
            continue
    f.close()
