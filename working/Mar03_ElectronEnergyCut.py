from copy import deepcopy
from CCD.working.Feb24_GammaReconAtExit import EnergyCosPhiParser

cutVars = ('cosPhi', 'energyExit', 'gammaEnergy', 'lineNumber')
print "10 um thick Nitro"
fn = "/Users/Yigong/GEANT4/work/BeamDiag_simple/build/out_GammaReconAtExit/nitro_10um.out"
Nitro_10um = EnergyCosPhiParser(fileName=fn, scatterWidth=10)
Nitro_10um.parse()
Nitro_10um.energyRecon()
Nitro_10um.addCutVar(cutVars)
Nitro_out = deepcopy(Nitro_10um)
Nitro_in = deepcopy(Nitro_10um)
keep = (Nitro_in.gammaEnergy > 1723.3)  & (Nitro_in.gammaEnergy < 1740.6)  
Nitro_in.cutData(keep)
Nitro_out.cutData(~keep)
ratioIn = []
ratioInOverTotal = []
total = len(Nitro_in.gammaEnergy) + len(Nitro_out.gammaEnergy)
total *= 1.
for eCut in range(0, 1600, 100):
    
    keepIn = Nitro_in.energyExit > eCut
    Nitro_in.cutData(keepIn)
    print 'keepIn has shape', keepIn.shape
    print 'Nitro_in has shape', Nitro_in.gammaEnergy.shape
    keepOut = Nitro_out.energyExit > eCut
    Nitro_out.cutData(keepOut)
    print 'keepOut has shape', keepOut.shape
    print 'Nitro_out has shape', Nitro_out.gammaEnergy.shape
    countsIn_tmp = len(Nitro_in.gammaEnergy) * 1.
    ratioInOverTotal.append(countsIn_tmp/(total))
    countsOut_tmp = len(Nitro_out.gammaEnergy) * 1.
    ratioIn.append(countsIn_tmp/(countsIn_tmp + countsOut_tmp))
    
fig, ax = plt.subplots(figsize=(10,6))
eCut = np.arange(0, 1600, 100)
ax.plot(eCut, ratioIn, 'b', label='Total: the number of events after e- energy cut')
ax.plot(eCut, ratioInOverTotal, 'b', label='Total: the number of events being detected')
ax.legend(loc=2)
ax.set_xlabel('Cut on electron energy (keV)')
ax.set_ylabel('Peak-to-Total ratio')
ax.set_ylim(ymin=0)
ax.grid()
ax.set_ylim((0,1))
fig.show()
    