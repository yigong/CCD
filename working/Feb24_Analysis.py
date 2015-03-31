from CCD.working.Feb24_GammaReconAtExit import EnergyCosPhiParser
import pickle

fig, ax = plt.subplots()
BINS = np.arange(0,4000,0.1)
cutVars = ('cosPhi', 'energyExit', 'gammaEnergy', 'lineNumber')
energyThreshold = 0

print "10 um thick Nitro"
fn = "/Users/Yigong/GEANT4/work/BeamDiag_simple/build/out_GammaReconAtExit/nitro_10um.out"
Nitro_10um = EnergyCosPhiParser(fileName=fn, scatterWidth=10)
Nitro_10um.parse()
Nitro_10um.energyRecon()
Nitro_10um.addCutVar(cutVars)
# # Nitro_10um = pickle.load(open('Nitro_10um.p', 'rb'))
Nitro_10um.cutData(Nitro_10um.energyExit > energyThreshold)
# Nitro_10um.hist(ax)
# Nitro_10um.GaussianFit(ax)
print 

print "50 um thick Nitro"
fn = "/Users/Yigong/GEANT4/work/BeamDiag_simple/build/out_GammaReconAtExit/nitro_50um.out"
Nitro_50um = EnergyCosPhiParser(fileName=fn, scatterWidth=50)
Nitro_50um.parse()
Nitro_50um.energyRecon()
Nitro_50um.addCutVar(cutVars)
#pickle.dump(Nitro_50um, open('Nitro_50um.p', 'wb'))
# Nitro_50um = pickle.load(open('Nitro_50um.p', 'wb'))
Nitro_50um.cutData(Nitro_50um.energyExit>energyThreshold)
# Nitro_50um.hist(ax, color='r')
# Nitro_50um.GaussianFit(ax, color='r')
print 
 
print "100 um thick Nitro"
fn = "/Users/Yigong/GEANT4/work/BeamDiag_simple/build/out_GammaReconAtExit/nitro_100um.out"
Nitro_100um = EnergyCosPhiParser(fileName=fn, scatterWidth=100)
Nitro_100um.parse()
Nitro_100um.energyRecon()
Nitro_100um.addCutVar(cutVars)
# #pickle.dump(Nitro_100um, open('Nitro_100um.p', 'wb'))
# # Nitro_100um = pickle.load(open('Nitro_100um.p', 'wb'))
Nitro_100um.cutData(Nitro_100um.energyExit>energyThreshold)
#figure()
Nitro_10um.hist(ax, color='b', bins=BINS, norm=False)
# Nitro_10um.GaussianFit(ax, color='b')
#figure()
# Nitro_50um.hist(ax, color='r', bins=BINS, norm=False)
# # Nitro_50um.GaussianFit(ax, color='r')
# #figure()
# Nitro_100um.hist(ax, color='y', bins=BINS, norm=False)
# Nitro_100um.GaussianFit(ax, color='y')
# for eThreshold in energyThresholds:
     



ax.set_xlabel('Gamma-ray Energy (keV)   binWidth = 0.1 keV')
ax.set_ylabel('Counts normalized to the maximum')
# ax.set_yscale('log')
# ax.set_ylim((1,1000))
#ax.set_title('electron energy > {:.0f} keV'.format(energyThreshold))
fig.show()
# fig.canvas.draw()

#