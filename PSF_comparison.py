from glob import glob
import re
from pickle import load, dump
import matplotlib.pyplot as plt
import numpy as np

fileList = glob('/Users/Yigong/Google Drive/Research/CCD/data/step_size/PSF_*[un]m_*P.p')
fig, ax = plt.subplots()
fig_alpha, ax_alpha = plt.subplots()

ll = dict()         # line style dict
ll['backP']  = '--'
ll['pixelP'] = '-'

lc = dict()
lc['10um'] = 'm'
lc['3um'] = 'c'
lc['1um'] = 'g'
lc['300nm'] = 'r'
lc['100nm'] = 'b'

lwd_x=1.5
lwd_a=1
for file in fileList:
    step = re.findall('\d+[un]m', file)[0]
    plane = re.findall('backP|pixelP', file)[0]
          

    [xWindow, alphaT, alphaM, betaM, xT, yT, xM, yM, ET, EM] = load(open(file, 'rb'))
#     if plane == 'backP':
#         H, xedges, yedges = np.histogram2d()
#     scatter(xT, xM)
#     xlabel('True x')
#     ylabel('Esti x')
#     title('%s_%s'%(step, plane))
#     scatter(yT, yM)
#     xlabel('True y')
#     ylabel('Esti y')
#     title('%s_%s'%(step, plane))
 
    # Eftr = ET > 1000.
#     xT = list(xT)
#     yT = list(yT)
#     xT_ftr = xT > -6500 and xT < 6500
#     yT_ftr = yT > 2000 and yT < 15000
#     ftr = xT_ftr and yT_ftr
    ftr = np.ones_like(EM)
    bins_xWindow = np.arange(-50, 50.01, 0.5)
    hist_xWindow, dummy = np.histogram(xWindow[ftr]/1e3, bins_xWindow)
    bins_alpha =  np.arange(-180, 540.01, 1.0)
    hist_alphaM, dummy = np.histogram(alphaM[ftr], bins_alpha)
    hist_alphaT, dummy = np.histogram(alphaT[ftr], bins_alpha)

    # plot
    ax.plot(bins_xWindow[1:], hist_xWindow, c=lc[step], lw=lwd_x, ls=ll[plane], label='%s_%s' %(step, plane))
    if (plane == 'backP') and (step == '100nm'):
        ax_alpha.plot(bins_alpha[1:], hist_alphaT, c=lc[step], label='alphaT', ls='-', lw=2)
    ax_alpha.plot(bins_alpha[1:], hist_alphaM, c=lc[step], lw=lwd_a, ls=ll[plane], label='%s_%s' %(step, plane))
    
# edit plot features

ax.legend()
ax.set_xlabel('Back-projected position at window (mm)\n bin_width = 0.5 mm')
ax.set_ylabel('Counts')
# ax.set_yscale('log')
ax.grid()

ax_alpha.legend()
ax_alpha.set_xlabel('Alpha (deg)')
# ax_alpha.set_yscale('log')
ax_alpha.grid()

plt.show()
