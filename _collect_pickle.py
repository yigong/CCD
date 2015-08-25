from glob import glob
import pickle
import numpy as np
from pylab import *
import traceback

try:
    pickle_directory = '/global/scratch/ygzhang/G4BeamDiag/track_recon/out/1200keV1um/pickle'
    alpha_linearReg = np.zeros(10001)
    alpha_median    = np.zeros_like(alpha_linearReg)
    for i in xrange(1, 10001):
        p_name = '%s/%d.p' %(pickle_directory, i)
        if glob(p_name):
            result = pickle.load(open(p_name, 'rb'))
            alpha_linearReg[i] = result['alpha_linearReg']
            alpha_median[i]    = result['alpha_median']
    bin_min = -150.
    bin_max = -90.
    bin_step = 1.
    bin_edges = np.arange(bin_min-0.5*bin_step, bin_max+0.51*bin_step, bin_step)
    bin_centers = np.arange(bin_min, bin_max+0.001, bin_step)
    alpha_linearReg_counts, _ = np.histogram(alpha_linearReg[1:], bins=bin_edges)
    alpha_median_counts, _    = np.histogram(alpha_median[1:],    bins=bin_edges)
    figure();
    plot(bin_centers, alpha_linearReg_counts, 'bs-', lw=2, label='linear regression')
    plot(bin_centers, alpha_median_counts, 'ro-', lw=2, label='median')
    legend()
    savefig('%s/LR_median.png' %pickle_directory)
except:
    traceback.print_exc()
