fig, (ax1, ax2) = plt.subplots(1,2)
ax1.hist(energy_list[12448:], bins=np.arange(0,1000), log=True)
ax1.set_xlabel('Energy (keV)')
ax2.hist(length_list[12448:], bins=np.arange(0,1000), log=True)
ax2.set_xlabel('Length (# pixels)') 
plt.show()