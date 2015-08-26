from scipy.io import loadmat
from glob import glob
import os
import pickle
import traceback
import numpy as np

try:
    matDir = '/global/scratch/ygzhang/G4BeamDiag/track_recon/out/1200keV1um/mat'
    os.chdir(matDir)
    alpha_BP = np.zeros(10001)
    for i in xrange(1, 10001):
        m_name = '%d.mat' % i
        if glob(m_name):
            a = loadmat(m_name, squeeze_me=True)
            alpha_BP[i] = a['alphaM']
    pickle.dump(alpha_BP, open('alpha.p', 'wb'))
except:
    traceback.print_exc()


