from CCD.working.G4Macro import G4Macro
from CCD.working.slurmGenerator import slurmGenerator
from CCD.working.G4Bash import G4Bash
from CCD.working.track_analysis import eData
from CCD.working.PositionSampling import PositionSampling

edata = eData('./recon_tracks_2hr.mat')
edata.betaMapping()
edata.alphaMapping()
edata.energyMapping()
p = PositionSampling(edata.alpha_mean_map, edata.beta_mean_map, edata.phi_mean_map, edata.energy_mean_map_1173, edata.energy_mean_map_1332)
p.sample()
for i in xrange(16):
    g4 = G4Macro('eSample_%s.mac'%(int(i)))
    g4.writeAll(p.energySample[i], p.betaSample[i])
    slm = slurmGenerator(i)
    slm.writeAll()

Macros = ['eSample_%s.mac'%(int(i)) for i in range(16)]
bash = G4Bash('eSample.sh')
bash.writeAll(Macros)
