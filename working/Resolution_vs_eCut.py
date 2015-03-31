from CCD.ElectronEnergyCut import ElectronEnergyCut

fig, ax = plt.subplots(figsize=(12,8))

fn = '/Users/Yigong/GEANT4/work/BeamDiag_simple/build/out_GammaReconAtExit/nitro_100um.out'

colors = ['b', 'r', 'y']
resolutions = [17.32, 86.6, 173.2]
labels = ['1% BW', '5% BW', '10% BW']

for (clr, resol, label) in zip(colors, resolutions, labels):
    Nitro10um = ElectronEnergyCut(fileName=fn, scatterWidth=100)
    Nitro10um.defineEnergyResolution(energyResolution=resol)
    Nitro10um.separateInOut()
    Nitro10um.cutData()
    Nitro10um.plot(ax, color=clr, labels=[label, ''])
ax.legend(loc='upper left')
ax.grid()
ax.set_xlabel('Electron energy cut (keV)')
ax.set_ylabel('Peak-to-total ratio')
ax.set_ylim((0,1))
plt.show()