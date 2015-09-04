from glob import glob
import sys
from PyBeamDiag.track_recon import ridge_follow
from PyBeamDiag.fit_func import gauss
import multiprocessing
from functools import partial
from scipy.io import loadmat
from scipy.optimize import curve_fit
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import pickle
import numpy as np

def main():
    args = sys.argv
    fits_dir = args[1]
    recon_dir = args[2]
    try:
        HTout = args[3]
    except IndexError:
        HTout = '/global/scratch/ygzhang/G4BeamDiag/track_recon/out/1200keV1um/alpha_px.mat'

    # recon LLSE
    f_files = glob('%s/*.fits' %(fits_dir))
    partial_func = partial(ridge_follow, outdir=recon_dir, plotflag=True, 
                           pickleflag=False)
    pool = multiprocessing.Pool(16)
    LLSE = pool.map(partial_func, f_files)
    pickle.dump(LLSE, open('%s/alpha_px.p' %(recon_dir), 'wb'))

    # load mat
    HT = loadmat(HTout, squeeze_me=True)['result_arr']
    HT[HT<90] = -90 - HT[HT<90] 
    HT[HT>90] = 270 - HT[HT>90]

    #plot
    binW = 1 # deg
    binCenters = np.arange(-150, -90.1, binW)
    binEdges = np.arange(binCenters[0]-0.5*binW, binCenters[-1]+0.51*binW, binW)
    LLSE_counts, _ = np.histogram(LLSE, binEdges)
    HT_counts, _   = np.histogram(HT, binEdges)
    fig = Figure()
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_subplot(111)
    ax.plot(binCenters, LLSE_counts, 'bo-', lw=2, label='LLSE')
    ax.plot(binCenters, HT_counts, 'rs-', lw=2, label='HT')

    # compute FWHM
    p0 = [10000, -120, 10]
    LLSE_fit, LLSE_cov = curve_fit(gauss, binCenters, LLSE_counts, p0)
    LLSE_fit_std = np.sqrt(np.diag(LLSE_cov))
    HT_fit, HT_cov = curve_fit(gauss, binCenters, HT_counts, p0)
    HT_fit_std = np.sqrt(np.diag(HT_cov))
    x = np.arange(-150, -90.01, 0.1)
    ax.plot(x, gauss(x, *(LLSE_fit)), 'b-', lw=1) 
    ax.plot(x, gauss(x, *(HT_fit)), 'r-', lw=1)

    
    # write result to file
    f = open('%s/Result.print' %(recon_dir), 'w')
    f.write('LLSE \n')
    f.write('    FWHM     = %.1f (%.2f)\n' %(LLSE_fit[-1]*2.35, LLSE_fit_std[-1]*2.35))
    f.write('    Centroid = %.1f (%.2f)\n' %(LLSE_fit[-2], LLSE_fit_std[-2]))
    f.write('    Total    = %s   \n' %(np.sum(LLSE_counts)))
    f.write('HT \n')
    f.write('    FWHM     = %.1f (%.2f)\n' %(HT_fit[-1]*2.35, HT_fit_std[-1]*2.35))
    f.write('    Centroid = %.1f (%.2f)\n' %(HT_fit[-2], HT_fit_std[-2]))
    f.write('    Total    = %s   \n' %(np.sum(HT_counts)))
    f.write('HTout: %s' %(HTout))
    f.close()

    # print figure
    canvas.print_figure('%s/LLSE_vs_HT.png' %(recon_dir))

if __name__ == '__main__':
    main()
