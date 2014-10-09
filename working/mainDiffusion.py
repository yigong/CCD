edata = eData('./recon_tracks_2hr.mat')
edata.cutData(edata.track_energy > 500)
edata.cutData(edata.ends_num == 2)
