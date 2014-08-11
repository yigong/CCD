import numpy as np
import matplotlib.pyplot as plt
from pylab import *

fig0, ax0 = plt.subplots()
idx_img = range(len(baseline))
ax0.plot(idx_img, baseline, color='blue')
ax0.set_ylabel(r"baseline[ADU]", fontsize=18, color='blue')
ax0_tw = ax0.twinx()
ax0_tw.plot(idx_img, noise, color='red')
ax0_tw.set_ylabel(r"noise[1$\sigma$]", fontsize=18, color='red')
