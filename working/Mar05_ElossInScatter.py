from CCD.working.Feb24_GammaReconAtExit import EnergyCosPhiParser
from matplotlib.colors import LogNorm
from matplotlib.cm import gray

cutVars = ('cosPhi', 'energyExit', 'gammaEnergy', 'lineNumber')

fn = '/Users/Yigong/GEANT4/work/BeamDiag_simple/build/out_GammaReconAtExit/nitro_10um.out'
nitro_10um =EnergyCosPhiParser(fileName=fn)
nitro_10um.parse()
nitro_10um.computeEnergyAngleDiff()
nitro_10um.energyRecon()
nitro_10um.addCutVar(cutVars)

fig, ax = plt.subplots()

energyBins = np.arange(0, 1700, 5)
diffBins = np.arange(0, 300, 1)
counts, energyEdges, diffEdges, img = ax.hist2d(nitro_10um.energyExit, nitro_10um.energyDiff, bins=[energyBins, diffBins], \
          normed=LogNorm, cmap=gray)


ax.scatter(nitro_10um.energyExit, nitro_10um.energyDiff, lw=0, alpha=0.5)
ax.set_xlabel('Electron energy when exiting the scatter (keV)')
ax.set_ylabel('Difference in energy (keV)')
ax.set_xlim((0, 1700))
ax.set_ylim((1e-2, 200))
ax.set_yscale('log')
fig.show()