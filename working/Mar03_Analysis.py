from CCD.working.Mar03_DopplerBroadening_Iodine import DopplerBroadening
fig, ax = plt.subplots()
bins = np.arange(-50, 50, 0.10)  

fn = "/Users/Yigong/GEANT4/work/DopplerBroadening/build/DB_Iodine_203keV.out"
I_203 = DopplerBroadening(fileName=fn)
I_203.parse()
I_203.sort()
I_203.hist(ax,bins)
ax.set_xlabel('$E_{GEANT4} - E_{Computed}$ (keV)     binWidth = 0.1 keV')
ax.set_ylabel('Counts')
ax.set_title('10^6 203 keV photons are incident on 1 mm Iodine. There are total of 4420 Compton events')


fig.show()