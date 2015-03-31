fig, ax = plt.subplots(figsize=(8,5))
BINS = np.arange(0, 4000, 1)

Nitro_10um.hist(ax, color='k', bins=BINS, norm=False)
# Nitro_10um.GaussianFit(ax, color='b')
# figure()
# Nitro_50um.hist(ax, color='r', bins=BINS, norm=False)
# Nitro_50um.GaussianFit(ax, color='r')
# figure()
# Nitro_100um.hist(ax, color='y', bins=BINS, norm=False)
# Nitro_100um.GaussianFit(ax, color='y')

#ax.set_xlabel('Gamma-ray Energy (keV)   binWidth = {} keV'.format(BINS[1]-BINS[0]))
#ax.set_ylabel('Counts')

# if zoom in
ax.set_xlim((1500, 2000))
# ticks = [item for item in np.arange(1710, 1751, 10)]
# labels = [str(item) for item in np.arange(1710, 1751, 10)]
# ax.set_xticks(ticks)
# ax.set_xticklabels(labels)
#ax.legend_.remove()
#fig.canvas.draw()
ax.legend_.remove()
fig.show()
