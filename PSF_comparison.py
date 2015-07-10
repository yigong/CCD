from glob import glob
import re
from pickle import load, dump
import matplotlib.pyplot as plt
import numpy as np

fileList = glob('/Users/Yigong/Google Drive/Research/CCD/data/track_recon_*[un]m.p')
fig, ax = plt.subplots()
for file in fileList:
    step = re.findall('\d+[un]m', file)[0]
    [x_window, x_CCD, y_CCD, alpha, beta] = load(open(file, 'rb'))
    bins = np.arange(-50, 50.01, 0.5)
    hist, dummy = np.histogram(x_window/1e3, bins)
    # plot
    ax.plot(bins[1:], hist, label='%s' %(step))
# edit plot features
ax.legend()
ax.set_xlabel('Back-projected position at window (mm)\n bin_width = 0.5 mm')
ax.set_ylabel('Counts')
plt.grid()
plt.show()
