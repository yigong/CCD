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

    # recon LLSE
    f_files = glob('%s/*.fits' %(fits_dir))
#     for f in f_files:
#         ridge_follow(f, outdir=recon_dir, plotflag=True, pickleflag=True)
    partial_func = partial(ridge_follow, outdir=recon_dir, plotflag=True, 
                           pickleflag=False)
    pool = multiprocessing.Pool(60)
    LLSE = pool.map(partial_func, f_files)
    # print LLSE
    pickle.dump(LLSE, open('%s/LLSE.p' %(recon_dir), 'wb'))
    
    # load mat
    HTout = '/global/scratch/ygzhang/G4BeamDiag/track_recon/out/1200keV1um/mat/alpha.p'
    HT = pickle.load(open('%s' %(HTout), 'rb'))
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
    HT_fit, HT_cov = curve_fit(gauss, binCenters, HT_counts, p0)
    x = np.arange(-150, -90.01, 0.1)
    ax.plot(x, gauss(x, *(LLSE_fit)), 'b-', lw=1) 
    ax.plot(x, gauss(x, *(HT_fit)), 'r-', lw=1)

    # write result to file
    f = open('%s/Result.print' %(recon_dir), 'w')
    f.write('LLSE \n')
    f.write('    FWHM = %.1f \n' %(LLSE_fit[-1]*2.35))
    f.write('    Total= %s \n' %(np.sum(LLSE_counts)))
    f.write('HT \n')
    f.write('    FWHM = %.1f \n' %(HT_fit[-1]*2.35))
    f.write('    Total= %s \n' %(np.sum(HT_counts)))
    f.close()

    # print figure
    canvas.print_figure('%s/LLSE_vs_HT.png' %(recon_dir))

if __name__ == '__main__':
    main()
