figEe, axEe = plt.subplots()
figEe.subplots_adjust(bottom=0.13)
axEe.hist(Ee, 100)
axEe.set_xlabel('$E_e$ (keV)', fontsize='xx-large')
axEe.set_ylabel('Counts')

figPhi, axPhi = plt.subplots()
figPhi.subplots_adjust(bottom=0.13)
axPhi.hist(PhiRadius/np.pi*180, 100)
axPhi.set_xlabel('$\phi$ (degree)', fontsize='xx-large')
axPhi.set_ylabel('Counts')

plt.show()