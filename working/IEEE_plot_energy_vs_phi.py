fig, ax = plt.subplots()
fig.subplots_adjust(bottom=0.15, left=0.14)
a = ax.scatter(edataAfterEnergy.phi)

binsEnergy = np.arange(0,2001,60)
a = ax.hist(edataAfterEnergy.phi, bins=binsEnergy, histtype='stepfilled',\
            align='left', color='b')
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Counts')
ax.set_xlim(0, 2000)
ax.set_ylim(0, 30)
fig.show()