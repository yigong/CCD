from CCD.working.Feb24_GammaReconAtExit import EnergyCosPhiParser

cutVars = ('cosPhi', 'energyExit', 'gammaEnergy', 'lineNumber')

fn = '/Users/Yigong/GEANT4/work/BeamDiag_simple/build/out_GammaReconAtExit/nitro_10um.out'
nitro_10um =EnergyCosPhiParser(fileName=fn)
nitro_10um.parse()
nitro_10um.computeEnergyAngleDiff()
nitro_10um.energyRecon()
nitro_10um.addCutVar(cutVars)

fig, ax = plt.subplots()
ax.scatter(np.arccos(nitro_10um.cosPhi)*180/np.pi, nitro_10um.angleDiff*180/np.pi, lw=0, alpha=0.5)

fig.show()