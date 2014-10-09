from CCD.working.track_analysis import eData
from CCD.working.PositionSampling import PositionSampling

edata = eData('./recon_tracks_2hr.mat')
edata.betaMapping()
edata.alphaMapping()
edata.energyMapping()
p = PositionSampling(edata.alpha_mean_map, edata.beta_mean_map, edata.phi_mean_map, edata.energy_mean_map_1173, edata.energy_mean_map_1332)
