fig, ax = plt.subplots(figsize=(8,6))
BINS = np.arange(1720,1745,0.1)
# 
# Nitro_10um.hist(ax, color='b', bins=BINS, norm=False)
# Nitro_50um.hist(ax, color='r', bins=BINS, norm=False)
# Nitro_100um.hist(ax, color='y', bins=BINS, norm=False)

ax.set_xlabel('Bin Width = {} keV'.format(0.1))
ax.set_ylabel('')
ax.set_xlim((1700, 1760))
ax.set_ylim((0, 80))
#ax.legend_.remove()
plt.show()