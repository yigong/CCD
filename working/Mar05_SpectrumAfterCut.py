from copy import deepcopy

print "10 um thick Nitro"
fn = "/Users/Yigong/GEANT4/work/BeamDiag_simple/build/out_GammaReconAtExit/nitro_10um.out"
Nitro = EnergyCosPhiParser(fileName=fn, scatterWidth=10)
Nitro.parse()
Nitro.energyRecon()
Nitro.addCutVar(cutVars)
fig, ax = plt.subplots(figsize=(12,8))
BINS = np.arange(0, 2500,2)
Nitro.hist(ax, color='r', bins=BINS, norm=False)

eCut = 800
keep = Nitro.energyExit > eCut
Nitro.cutData(keep)
Nitro.hist(ax, color='b', bins=BINS, norm=False)
#ax.set_yscale('log')
ax.legend(['Before the 800 keV cut', 'After the cut'], loc=2)
ax.set_xlabel('Gamma-ray Energy (keV)   binWidth = {} keV'.format(BINS[1]-BINS[0]))
ax.set_ylabel('Counts [log scale]')
ax.set_xlim((1500,2000))
ax.set_ylim((0,100))
fig.show()

